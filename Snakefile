# VariantHunter Version 0.5.0

configfile: "config.yaml"

import os

# open sample name file and read the content in a list
with open(config["sample_name_path"], 'r') as file:
    SAMPLE_NAMES = [sample_name.rstrip() for sample_name in file.readlines()]

rule all:
    input:
        "output/5_results/cluster_lookup.txt"

## Reference ##################################################################

# Remove illegal characters from reference names and index the reference
rule prepare_index_reference:
    input:
        config["ref_path"]
    output:
        fasta="output/1_reference/ref.fasta",
        c="output/1_reference/ref.comp.b",
        fai="output/1_reference/ref.fasta.fai",
        name="output/1_reference/ref.name",
        old_names=temp("output/1_reference/ref.old_name.txt"),
        converstion_file="output/1_reference/gene_name_conversion.tsv"
    shell:
        '''
        cp {input} {output.fasta};
        sed -i -e "s/[;|:|,| |(|)|'|+]/_/g" -e "s/-/_/g" -e "s/\//_/g" -e "s/\./_/g"  -e "s/|/_/g"  {output.fasta};
        sed -i 's/__/_/g'  {output.fasta};
        sed -i 's/__/_/g'  {output.fasta};

        kma index -i {output.fasta} -o output/1_reference/ref;
        samtools faidx {output.fasta};
        grep -e ">" {input} | sed 's/>//g'  > {output.old_names};
        paste {output.old_names} {output.name} > {output.converstion_file}

        '''

## Mapping #########################################################################

# Map samples to resfinder database, to create consensus fasta and other files
rule mapping:
    input:
        ref="output/1_reference/ref.comp.b"
    output:
        bam_raw=temp("output/2_mapping/bam_raw/{sample_name}.bam"),
        res_file="output/2_mapping/{sample_name}.res",
        bam_sort="output/2_mapping/{sample_name}.bam",
        bai="output/2_mapping/{sample_name}.bam.bai"
    params:
        S=config["seq_path"] + "{sample_name}" + config["singleton_suf"],
        out_dir="output/2_mapping/{sample_name}",
        PE=config["seq_path"] + "{sample_name}" + config["for_rev_suf"],
        ref="output/1_reference/ref",
        conclave=config["conclave_verson"],
        ali_settings=config["alignment_setting_for_kma"]

    shell:
        "kma -i {params.S} -ipe {params.PE} -o {params.out_dir} -t_db {params.ref} "
        " -1t1 -sam 2096 -ConClave {params.conclave} {params.ali_settings} -nc -nf | samtools view -bS > {output.bam_raw};"
        "samtools sort -o {output.bam_sort} {output.bam_raw};"
        "samtools index {output.bam_sort};"


## vcf ################################################################################

rule bcftools_mpileup_call:
    input:
        bam="output/2_mapping/{sample_name}.bam",
        ref="output/1_reference/ref.fasta"
    output:
        vcf="output/3_vcf/1_raw/{sample_name}.raw.vcf.gz"
    params:
        max_dp="--max-depth 1000000", # At a position, read maximally INT reads per input file. 
        MQ=config["minimum_map_Q"],   # Minimum mapping quality for an alignment to be used. 
        BQ=config["minimum_base_Q"],  # Minimum base quality for a base to be considered.
        p=config["p_value"],          # Ignore variant if the frequency of the ref allele  >= p
        ploidy="1",                   # 1 means haploid.
        orphan=config["orphan_reads"] # Choose whether or not to include orphan reads

    shell:
        "bcftools mpileup -Ou -f {input.ref} {params.max_dp} --min-MQ {params.MQ} --min-BQ {params.BQ} {params.orphan} "
        " -a INFO/AD  {input.bam} | "
        "bcftools call -c -p {params.p} --ploidy {params.ploidy} -Oz -o {output.vcf};"
        

# Normalize and left-align indels and remove dublicate entries of the same indel.
rule bcftools_norm_vcf:
    input:
        vcf_raw="output/3_vcf/1_raw/{sample_name}.raw.vcf.gz",
        ref="output/1_reference/ref.fasta"
    output:
        norm_vcf="output/3_vcf/2_norm/{sample_name}.norm.vcf.gz"
    shell:
        '''
        bcftools norm --rm-dup indels -f {input.ref} {input.vcf_raw} -Oz -o {output.norm_vcf}
        '''
        
rule bcftools_filter_vcf:
    input:
        vn="output/3_vcf/2_norm/{sample_name}.norm.vcf.gz"
    output:
        vf="output/3_vcf/3_filter_norm/{sample_name}.filter.norm.vcf.gz",
        vf_index="output/3_vcf/3_filter_norm/{sample_name}.filter.norm.vcf.gz.csi"
    params:
        min_AD=config["minimum_ALT_AD"],   # minimum ALT allele depth to be accepted.
        min_AF=config["minimum_ALT_AF"],    # minimum ALT allele freq to be accepted.
        min_strand=config["min_ALT_strand_prop_F_and_R"]
    shell:
        '''
        bcftools filter -i '(ALT="." || (DP4[2] >= AD[1]*{params.min_strand} && \
        DP4[3]  >= AD[1]*{params.min_strand} && AD[1] = max(AD) && AD[1] >= {params.min_AD} && \
        AD[1] >= (INFO/DP*{params.min_AF} )))' {input.vn} -Oz -o {output.vf};
        bcftools index -o {output.vf_index} {output.vf}
        '''

## Sample consensus sequence wrangling  ##########################################

rule bcftools_consensus_raw:
    input:
        vf="output/3_vcf/3_filter_norm/{sample_name}.filter.norm.vcf.gz",
        vf_index="output/3_vcf/3_filter_norm/{sample_name}.filter.norm.vcf.gz.csi",
        ref="output/1_reference/ref.fasta"
    output:
        raw_consensus="output/4_induvidual_sample_seqs/1_raw/{sample_name}.fasta"
    shell:
        '''
        bcftools consensus --absent N  -f {input.ref} {input.vf} > {output.raw_consensus}
        '''


# Make coverage statistics
rule samtools_coverage:
    input:
        bam_sort="output/2_mapping/{sample_name}.bam"
    output:
        cov="output/4_induvidual_sample_seqs/0_coverage_stats/{sample_name}.cov"
    params:
        min_BQ=config["minimum_base_Q"]
    shell:
        "samtools coverage --no-header --min-BQ {params.min_BQ}  {input.bam_sort} > {output.cov}"



# Filter  consensus sequences on coverage
rule filter_coverage:
    input:
        r_c="output/4_induvidual_sample_seqs/1_raw/{sample_name}.fasta",
        cov="output/4_induvidual_sample_seqs/0_coverage_stats/{sample_name}.cov"
    output:
        hc_c="output/4_induvidual_sample_seqs/2_high_cov/{sample_name}.fasta",
        hc_list=temp("output/4_induvidual_sample_seqs/2_high_cov/{sample_name}.hc_list")
    params:
        min_cov=config["min_template_cov"],    # minimum template coverage in percent to accept consensus
        min_meandepth=config["min_meandepth"]
    run:
        shell("awk '$6 >= {params.min_cov} && $7 >= {params.min_meandepth}' {input.cov} | awk '{{print $1}}'  > {output.hc_list}")
        if os.path.getsize(output[1]) > 0:  # Test that high coverage list file is not empty.
            shell("seqkit grep -n -f {output.hc_list} {input.r_c} > {output.hc_c}")
        else:
            shell("touch {output.hc_c}")    # Create empty fasta file if high coverage list is empty
        

# Filter on percentage of Ns in sequences
rule seq_cleaner_filter_Ns:
    input:
        "output/4_induvidual_sample_seqs/2_high_cov/{sample_name}.fasta"
    output:
        "output/4_induvidual_sample_seqs/3_high_cov_n/{sample_name}.fasta"
    params:
        l="0",  # minimum length of sequences. Set to 0 since coverage is filtered upstreams in the pipeline.
        p=config["percentage_unknown_nuc_allowed"]   # max percentage Ns allowed in sequences
    shell:
        "cd output/4_induvidual_sample_seqs/2_high_cov;"
        "python3 ../../../scripts/sequence_cleaner.py  {wildcards.sample_name}.fasta  {params.l} {params.p};"
        "mv clear_{wildcards.sample_name}.fasta ../3_high_cov_n/{wildcards.sample_name}.fasta"

rule make_list_of_templates_used:
    input:
        expand("output/4_induvidual_sample_seqs/3_high_cov_n/{sample_name}.fasta", sample_name=SAMPLE_NAMES)
    output:
        temp("output/5_results/template_names.txt")
    shell:
        "cat output/4_induvidual_sample_seqs/3_high_cov_n/*.fasta | grep '>' | sed 's/>//g' | sort -u > {output}"


rule perl_sample_id_to_header_sample_fasta:
    input: 
        "output/4_induvidual_sample_seqs/3_high_cov_n/{sample_name}.fasta"
    output:
        "output/4_induvidual_sample_seqs/4_high_cov_n_ID/{sample_name}.fasta"
    params:
         sample_id= lambda wildcards: wildcards.sample_name
    shell:
        '''   
        perl -p -e 's/^>/>{params.sample_id}+/g' {input}  > {output}
        '''     
## Pool Samples ############################################################################

rule cat_variant_fastas:
    input:
        gene_fasta=expand("output/4_induvidual_sample_seqs/4_high_cov_n_ID/{sample_name}.fasta", sample_name=SAMPLE_NAMES)
    output:
        "output/5_results/all_consensus_sequences.fasta"
    shell:
        '''
        cat output/4_induvidual_sample_seqs/4_high_cov_n_ID/*.fasta   >   {output}
        '''    

rule add_ref_sequences:
    input:
        ref_names="output/5_results/template_names.txt",
        cons_fasta="output/5_results/all_consensus_sequences.fasta",
        ref_fasta="output/1_reference/ref.fasta"
    output:
        both_fasta="output/5_results/all_consensus_and_ref_sequences.fasta",
        ref_subset_fasta=temp("output/5_results/ref_sequences.fasta")
    shell:
        "seqkit grep -n -f {input.ref_names} {input.ref_fasta} | perl -p -e 's/^>/>Ref+/g' > {output.ref_subset_fasta};"
        "cat {input.cons_fasta} {output.ref_subset_fasta}  >   {output.both_fasta}"


# Merge dublicate sequences so that all sequences are unique. Also merges the headers
rule DupRemover:
    input:
        "output/5_results/all_consensus_and_ref_sequences.fasta"
    output:
        temp("output/5_results/all_consensus_and_ref_sequences_dublicants_merged_ML.fasta")
    shell:
        "python3 scripts/DupRemover.py -i {input} -o {output};"


rule single_line_fastas:
    input:
        "output/5_results/all_consensus_and_ref_sequences_dublicants_merged_ML.fasta"
    output:
        temp("output/5_results/all_consensus_and_ref_sequences_dublicants_merged.fasta")
    shell:
        '''
        cat {input} | awk '/^>/ {{printf("\\n%s\\n",$0);next; }} {{ printf("%s",$0);}}  END {{printf("\\n");}}' | tail -n+2  > {output}
        '''


# Give the sequences proper names and create 2 metadata files
rule python3_create_variant_ID_and_metafiles:
    input:
        fasta="output/5_results/all_consensus_and_ref_sequences_dublicants_merged.fasta",
        gene_name_conversion_file="output/1_reference/gene_name_conversion.tsv"
    output:
        fasta_id="output/5_results/all_consensus_and_ref_sequences_dublicants_merged_ID.fasta",
        raw_tsv=temp("output/5_results/raw_sequence_sample_metadata.tsv"),
        tsv="output/5_results/sequence_sample_metadata.tsv",
        raw_template_metadata_file =temp("output/5_results/raw_template_metadata.tsv"),
        template_metadata_file ="output/5_results/template_metadata.tsv"
    run:
    ## Create ID fasta and sequence sample metadata
        variant_meta = open("output/5_results/raw_sequence_sample_metadata.tsv", 'w')
        fasta_id_file = open("output/5_results/all_consensus_and_ref_sequences_dublicants_merged_ID.fasta",'w')
        fasta_file = open("output/5_results/all_consensus_and_ref_sequences_dublicants_merged.fasta","r")
        variant_count = {}
        refs_found=set()
        for line in fasta_file:
            if line.startswith(">"):
                line = line.strip(">")
                last_entry = line.split("|")[-1]
                gene_name= last_entry.split("+")[1].strip()  
                if last_entry.startswith("Ref"):
                    varname = gene_name + ".R"
                    if len(line.split("|")) > 1:        # if the reference was found in any actual samples and not just added.
                        refs_found.add(gene_name)
                else:
                    if gene_name not in variant_count:
                        variant_count[gene_name] = 1
                    else: 
                        variant_count[gene_name] += 1
                    varname = gene_name + ".v" + str(variant_count[gene_name])                    
                new_header= ">" + varname
                print(new_header, file =fasta_id_file)
                sample_and_gene_list = line.split("|") # split the string into the different headers
                sample_list = [i.split("+", 1)[0] for i in sample_and_gene_list] # make list of sample IDs
                if sample_list == ["Ref"]:
                    sample_list = ["*"]
                if "Ref" in sample_list:
                    sample_list.remove("Ref")
                if sample_list == ["*"]:
                    length = "0"
                else:            
                    length = str(len(sample_list))    
            else:
                print(line.strip(), file =fasta_id_file)
                metafile_line = varname + "\t" + length + "\t" + line.strip() +  "\t" + ','.join(sample_list)
                print(metafile_line, file =variant_meta)
        variant_meta.close()        
        fasta_id_file.close()
        fasta_file.close()
        shell("sort -k 2 -n -r  {output.raw_tsv} > {output.tsv}")  # sort by number of samples the sequence was found in
        shell("sed -i '1s/^/#Sequence_name\tNumber_of_samples\tSequence \tSample_list\\n/' {output.tsv}")           # Add Header to meta file
    ## Create Template metadata file
        gene_name_conversion_file = open("output/1_reference/gene_name_conversion.tsv","r")
        raw_template_metadata_file = open("output/5_results/raw_template_metadata.tsv","w")
        for line in gene_name_conversion_file:
            line = line.strip()
            new_name = line.split("\t")[1]
            if new_name in refs_found:
                ref_found = "yes"
                if new_name in variant_count:
                    found_vars=str(variant_count[new_name] + 1)
                else:
                    found_vars="1"
            else:
                ref_found = "no"
                if new_name in variant_count:
                    found_vars=str(variant_count[new_name])
                else:
                    found_vars="0"
            new_line= line + "\t" + found_vars + "\t" + ref_found
            print(new_line, file = raw_template_metadata_file)
        gene_name_conversion_file.close()
        raw_template_metadata_file.close()
        shell("sort -k 3 -n -r  {output.raw_template_metadata_file} > {output.template_metadata_file}")
        shell("sed -i '1s/^/#Official_gene_name\tSoftware_friendly_gene_name\tNumber_of_variants_found\tExact_Reference_sequence_found?\\n/' {output.template_metadata_file}") 

## Clusters ###############################################################################

rule cd_hit:
    input:
        "output/5_results/all_consensus_and_ref_sequences_dublicants_merged_ID.fasta"
    output:
        clus="output/5_results/cd_hit.clstr",
        genes=temp("output/5_results/cd_hit")
    params:
        i= config["cluster_identity_cd_hit"],
        w= config["cluster_wordsize_cd_hit"]
    shell:
        '''
        cd-hit-est -i {input} -o {output.genes} -M 0 -d 0 -c {params.i} -n {params.w} -d 0 -T 8 -s 0.5 -sc 1 -g 1 -mask N
        '''

rule make_cd_hit_fastas:
    input:
        clus="output/5_results/cd_hit.clstr",
        seq="output/5_results/all_consensus_and_ref_sequences_dublicants_merged_ID.fasta"
    output:
        temp(dynamic("output/5_results/cd-hit_fastas/{cluster_id}"))
    params:
        fasta_outdir="output/5_results/cd-hit_fastas",
        min_cluster_size = config["minimum_cluster_size"]
    shell:
        "perl scripts/make_multi_seq.pl {input.seq} {input.clus} {params.fasta_outdir} {params.min_cluster_size}"


# Some new variants may have ended up in an other cluster than their reference sequence.
# Add the reference sequence to these clusters before alignments.

rule python3_seqkit_add_refs_to_clusters_that_lost_refs:
    input:
        c_fasta_raw="output/5_results/cd-hit_fastas/{cluster_id}",
        all_seq_fasta="output/5_results/all_consensus_and_ref_sequences_dublicants_merged_ID.fasta"
    output:
        list=temp("output/5_results/ref_add_lists/{cluster_id}.ref_add.txt"),
        complete_c_fasta="output/5_results/cluster_fastas/{cluster_id}.fasta"
    run:
        fasta_in = open(input[0],"r")
        ref_set = set()
        var_set = set()
        for line in fasta_in:
            if line.startswith(">"):    # if it is a fasta header
                row = line[1:].strip()   # remove the >
                gene_name = row.split(".")[0] # extract gene name
                if row.split(".")[1] is "R":  # if it is a ref sequence
                    ref_set.add(gene_name)     # save the gene name in ref set
                else:                          # if it is not a ref sequence
                    var_set.add(gene_name)     # save the gene name in variant set              
        fasta_in.close()
        missing_ref = var_set.difference(ref_set)  # find the var genes that does not have their ref sequence in the cluster
        list_out= open(output[0],"w")
        for item in missing_ref:
            print(item + ".R", file=list_out)            
        list_out.close()
        shell("cp {input.c_fasta_raw}  {output.complete_c_fasta}")
        if len(missing_ref) > 0:  # if there are missing ref fastas, then add the missing ones to the input fasta
            shell("seqkit grep --line-width 0 -n -f {output.list} {input.all_seq_fasta} >> {output.complete_c_fasta}") 


rule align:
    input:
        fasta="output/5_results/cluster_fastas/{cluster_id}.fasta"
    output:
        fasta="output/5_results/cluster_alignments/{cluster_id}.fasta"
    shell:
        "mafft --auto {input.fasta} > {output.fasta}"


rule fasttree:
    input:
        "output/5_results/cluster_alignments/{cluster_id}.fasta"
    output:
        "output/5_results/cluster_trees/{cluster_id}.tree"
    shell:
        "FastTree -gtr -nt  {input} > {output}"

rule snp_dists:
    input:
        "output/5_results/cluster_alignments/{cluster_id}.fasta"
    output:
        "output/5_results/cluster_snp_dist/{cluster_id}.tsv"
    params:
        options =config["snp_dists_options"] 
    shell:
        "snp-dists {params.options} {input} > {output}"


rule make_cluster_sequence_sample_metadata:
    input:
        tree="output/5_results/cluster_trees/{cluster_id}.tree",
        tsv="output/5_results/sequence_sample_metadata.tsv",
        dist="output/5_results/cluster_snp_dist/{cluster_id}.tsv",
        fasta="output/5_results/cluster_fastas/{cluster_id}.fasta"
    output:
        tsv="output/5_results/cluster_meta/{cluster_id}.tsv",
        cc=temp("output/5_results/cluster_content/{cluster_id}.txt"),
        list=temp("output/5_results/cluster_list/{cluster_id}.txt")
    shell:
        "grep '>' {input.fasta} | sed 's/>//g' > {output.list};"
        "grep -wFf {output.list}  {input.tsv} > {output.tsv};"
        "grep '\.R' {output.list} | sed 's/\.R//g' | sed 's/^/{wildcards.cluster_id}    /' > {output.cc} "


rule make_cluster_lookup_file:
    input:
        dynamic("output/5_results/cluster_content/{cluster_id}.txt")
    output:
        "output/5_results/cluster_lookup.txt"
    shell:
        "cat output/5_results/cluster_content/*.txt | sort -n > {output}"
