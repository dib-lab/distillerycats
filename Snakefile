import pandas as pd
import feather
from sourmash import signature
import glob
import os
from collections import Counter

m = pd.read_csv("inputs/test_metadata.csv", header = 0)
SAMPLES = m['sample'].unique().tolist()

rule all:
    input:
        "outputs/comp/all_filt_comp.csv",
        "outputs/hash_tables/normalized_abund_hashes_wide.feather",
        #"outputs/hash_tables/all_unnormalized_abund_hashes_wide.feather",
        #"outputs/gather/vita_vars.csv",
        #"aggregated_checkpoints/aggregate_spacegraphcats_gather_matches.txt",
        #"aggregated_checkpoints/aggregate_spacegraphcats_gather_matches_plass.txt"

########################################
## PREPROCESSING
########################################

rule download_fastq_files:
    output: "inputs/raw/{sample}.fastq.gz",
    run:
        row = m.loc[m['sample'] == wildcards.sample]
        fastq = row['download'].values
        fastq = fastq[0]
        shell("wget -O {output} {fastq}")


rule download_human_db:
    output: "inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz"
    shell:'''
    wget -O {output} https://osf.io/84d59/download
    '''

rule remove_host:
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
    output:
        r = 'outputs/bbduk/{sample}.nohost.fq.gz',
        human='outputs/bbduk/{sample}.human.fq.gz',
    input: 
        r = 'inputs/raw/{sample}.fastq.gz',
        human='inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    conda: 'bbmap.yml'
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r} out={output.r} outm={output.human} k=31 ref={input.human}
    '''

rule kmer_trim_reads:
    input: 'outputs/bbduk/{sample}.nohost.fq.gz'
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    params: mem="20e9"
    conda: 'sourmash.yml'
    shell:'''
    trim-low-abund.py --gzip -C 3 -Z 18 -M {params.mem} -V {input} -o {output}
    '''

rule compute_signatures:
    input: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs/sigs/{sample}.sig"
    conda: 'sourmash.yml'
    shell:'''
    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o {output} {input}
    '''

########################################
## Filtering and formatting signatures
########################################

rule get_greater_than_1_filt_sigs:
    input: expand("outputs/sigs/{sample}.sig", sample = SAMPLES) 
    output: "outputs/filt_sig_hashes/greater_than_one_count_hashes.txt"
    run:
        # Determine the number of hashes, the number of unique hashes, and the number of
        # hashes that occur once across 954 IBD/control gut metagenomes (excludes the 
        # iHMP). Calculated for a scaled of 2k. 9 million hashes is the current 
        # approximate upper limit with which to build a sample vs hash abundance table 
        # using my current methods.

        files = input

        all_mins = []
        for file in files:
            if os.path.getsize(file) > 0:
                sigfp = open(file, 'rt')
                siglist = list(signature.load_signatures(sigfp))
                loaded_sig = siglist[1]
                mins = loaded_sig.minhash.get_mins() # Get the minhashes 
                all_mins += mins

        counts = Counter(all_mins) # tally the number of hashes

        # remove hashes that occur only once
        for hashes, cnts in counts.copy().items():
            if cnts < 2:
                counts.pop(hashes)

        # write out hashes to a text file
        with open(str(output), "w") as f:
            for key in counts:
                print(key, file=f)


rule convert_greater_than_1_hashes_to_sig:
    input: "outputs/filt_sig_hashes/greater_than_one_count_hashes.txt"
    output: "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig"
    conda: 'sourmash.yml'
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name greater_than_one_count_hashes --filename {input} {input}
    '''

rule filter_signatures_to_greater_than_1_hashes:
    input:
        filt_sig = "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig",
        sigs = "outputs/sigs/{sample}.sig"
    output: "outputs/filt_sigs/{sample}_filt.sig"
    conda: 'sourmash.yml'
    shell:'''
    sourmash sig intersect -o {output} -A {input.sigs} -k 31 {input.sigs} {input.filt_sig}
    '''

rule name_filtered_sigs:
    input: "outputs/filt_sigs/{sample}_filt.sig"
    output: "outputs/filt_sigs_named/{sample}_filt_named.sig"
    conda: 'sourmash.yml'
    shell:'''
    sourmash sig rename -o {output} -k 31 {input} {wildcards.sample}_filt
    '''

rule convert_greater_than_1_signatures_to_csv:
    input: "outputs/filt_sigs_named/{sample}_filt_named.sig"
    output: "outputs/filt_sigs_named_csv/{sample}_filt_named.csv"
    conda: 'sourmash.yml'
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_abund_table_long_normalized:
    input: 
        expand("outputs/filt_sigs_named_csv/{sample}_filt_named.csv", sample = SAMPLES)
    output: csv = "outputs/hash_tables/normalized_abund_hashes_long.csv"
    conda: 'r.yml'
    script: "scripts/normalized_hash_abund_long.R"

rule make_hash_abund_table_wide:
    input: "outputs/hash_tables/normalized_abund_hashes_long.csv"
    output: "outputs/hash_tables/normalized_abund_hashes_wide.feather"
    run:
        import pandas as pd
        import feather
        
        ibd = pd.read_csv(str(input), dtype = {"minhash" : "int64", "abund" : "float64", "sample" : "object"})
        ibd_wide=ibd.pivot(index='sample', columns='minhash', values='abund')
        ibd_wide = ibd_wide.fillna(0)
        ibd_wide['sample'] = ibd_wide.index
        ibd_wide = ibd_wide.reset_index(drop=True)
        ibd_wide.columns = ibd_wide.columns.astype(str)
        ibd_wide.to_feather(str(output)) 


########################################
## Random forests & optimization
########################################

rule install_pomona:
    input: "outputs/hash_tables/normalized_abund_hashes_wide.feather"
    output:
        pomona = "outputs/vita_rf/pomona_install.txt"
    conda: 'rf.yml'
    script: "scripts/install_pomona.R"

rule vita_var_sel_rf:
    input:
        info = "inputs/working_metadata.tsv", 
        feather = "outputs/hash_tables/normalized_abund_hashes_wide.feather",
        pomona = "outputs/vita_rf/pomona_install.txt"
    output:
        vita_rf = "outputs/vita_rf/vita_rf.RDS",
        vita_vars = "outputs/vita_rf/vita_vars.txt",
        ibd_novalidation = "outputs/vita_rf/ibd_novalidation_filt.csv",
        ibd_novalidation_diagnosis = "outputs/vita_rf/bd_novalidation_filt_diagnosis.txt",
        ibd_validation = "outputs/vita_rf/ibd_validation_filt.csv"
    conda: 'rf.yml'
    script: "scripts/vita_rf.R"

rule tune_rf:
    input:
        ibd_novalidation = "outputs/vita_rf/ibd_novalidation_filt.csv",
        ibd_novalidation_diagnosis = "outputs/vita_rf/bd_novalidation_filt_diagnosis.txt" 
    output:
        optimal_rf = "outputs/optimal_rf/optimal_ranger.RDS",
        pred_test = "outputs/optimal_rf/pred_test_tab.txt",
        pred_train = "outputs/optimal_rf/pred_train_tab.txt"
    conda: 'rf.yml'
    script: "scripts/tune_rf.R"

rule validate_rf:
    input:
        optimal_rf = "outputs/optimal_rf/optimal_ranger.RDS",
        ibd_validation = "outputs/vita_rf/ibd_validation_filt.csv",
        info = "inputs/working_metadata.tsv" 
    output:
        pred_srp = "outputs/rf_validation/pred_srp057027.txt",
        pred_srp_df = "outputs/rf_validation/pred_srp057027.csv",
        pred_prjna = "outputs/rf_validation/pred_prjna385949.txt",
        pred_prjna_df = "outputs/rf_validation/pred_prjna385949.csv"
    conda: "rf.yml"
    script: "scripts/validate_rf.R"



########################################
## Predictive hash characterization
########################################

rule convert_vita_vars_to_sig:
    input: "outputs/vita_rf/vita_vars.txt"
    output: "outputs/vita_rf/vita_vars.sig"
    conda: "sourmash.yml"
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name vita_vars --filename {input} {input}
    '''

rule download_gather_almeida:
    output: "inputs/gather_databases/almeida-mags-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/5jyzr/download
    '''

rule untar_almeida:
    output: "inputs/gather_databases/almeida-mags-k31.sbt.json"
    input: "inputs/gather_databases/almeida-mags-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_pasolli:
    output: "inputs/gather_databases/pasolli-mags-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/3vebw/download
    '''

rule untar_pasolli:
    output: "inputs/gather_databases/pasolli-mags-k31.sbt.json"
    input: "inputs/gather_databases/pasolli-mags-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_nayfach:
    output: "inputs/gather_databases/nayfach-k31.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/y3vwb/download
    '''

rule untar_nayfach:
    output: "inputs/gather_databases/nayfach-k31.sbt.json"
    input: "inputs/gather_databases/nayfach-k31.tar.gz"
    params: outdir="inputs/gather_databases"
    shell:'''
    tar xf {input} -C {params.outdir}
    '''

rule download_gather_genbank:
    output: "inputs/gather_databases/genbank-d2-k31.tar.gz"
    shell:'''
    wget -O {output} https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/genbank-d2-k31.tar.gz
    '''

rule untar_genbank:
    output: "inputs/gather_databases/genbank-d2-k31.sbt.json"
    input:  "inputs/gather_databases/genbank-d2-k31.tar.gz"
    params: outdir = "inputs/gather_databases"
    shell: '''
    tar xf {input} -C {params.outdir}
    '''

rule gather_vita_vars:
    input:
        sig="outputs/vita_rf/vita_vars.sig",
        db1="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db3="inputs/gather_databases/nayfach-k31.sbt.json",
        db4="inputs/gather_databases/pasolli-mags-k31.sbt.json"
    output: 
        csv="outputs/gather/vita_vars.csv",
        matches="outputs/gather/vita_vars.matches",
        un="outputs/gather/vita_vars.un"
    conda: 'sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db1} {input.db4} {input.db3} {input.db2}
    '''

rule download_gather_match_genomes:
    output: "outputs/gather/gather_genomes.tar.gz"
    shell:'''
    wget -O {output} https://osf.io/ungza/download
    '''

rule untar_gather_match_genomes:
    output:  directory("outputs/gather_genomes/")
    input:"outputs/gather/gather_genomes.tar.gz"
    params: outdir = "outputs/gather_genomes"
    shell:'''
    mkdir -p {params.outdir}
    tar xf {input} -C {params.outdir}
    '''

rule gtdbtk_gather_matches:
    """
    this rule require the gtdbtk databases. The tool finds the database by 
    using a path specified in a file in the environment. I predownloaded the 
    databases and placed them in the required location.
    The path is in this file:
    .snakemake/conda/9de8946b/etc/conda/activate.d/gtdbtk.sh
    """
    input: directory("outputs/gather_genomes/")
    output: "outputs/gtdbtk/gtdbtk.bac120.summary.tsv"
    params:  outdir = "outputs/gtdbtk"
    conda: "gtdbtk.yml"
    shell:'''
    gtdbtk classify_wf --genome_dir {input} --out_dir {params.outdir} --cpus 8 
    '''

checkpoint spacegraphcats_gather_matches:
    input: 
        query = directory("outputs/gather_genomes/"),
        conf = expand("inputs/sgc_conf/{sample}_r1_conf.yml", sample = SAMPLES),
        reads = expand("outputs/abundtrim/{sample}.abundtrim.fq.gz", sample = SAMPLES)
    output: 
        directory(expand("outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/", sample = SAMPLES))
        #"outputs/sgc_genome_queries/{library}_k31_r1_search_oh0/{gather_genome}.cdbg_ids.reads.fa.gz",
    params: outdir = "outputs/sgc_genome_queries"
    conda: "spacegraphcats.yml"
    shell:'''
    python -m spacegraphcats {input.conf} extract_contigs extract_reads --nolock --outdir={params.outdir}  
    '''

rule calc_sig_nbhd_reads:
    input: "outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{gather_genome}.fna.cdbg_ids.reads.fa.gz"
    output: "outputs/nbhd_read_sigs/{sample}/{gather_genome}.cdbg_ids.reads.sig"
    conda: "sourmash.yml"
    shell:'''
    sourmash compute -k 21,31,51 --scaled 2000 --track-abundance -o {output} --merge {wildcards.sample}_{wildcards.gather_genome} {input}
    '''

def aggregate_spacegraphcats_gather_matches(wildcards):
    # checkpoint_output produces the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.spacegraphcats_gather_matches.get(**wildcards).output[0]    
    file_names = expand("outputs/nbhd_read_sigs/{sample}/{gather_genome}.cdbg_ids.reads.sig",
                        sample = SAMPLES, 
                        gather_genome = glob_wildcards(os.path.join(checkpoint_output, "{gather_genome}.fna.cdbg_ids.reads.fa.gz")).gather_genome)
    return file_names


rule aggregate_signatures:
    input: aggregate_spacegraphcats_gather_matches
    output: "aggregated_checkpoints/aggregate_spacegraphcats_gather_matches.txt"
    shell:'''
    touch {output}
    '''

rule plass_nbhd_reads:
    input: "outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{gather_genome}.fna.cdbg_ids.reads.fa.gz"
    output: "outputs/nbhd_read_plass/{sample}/{gather_genome}.cdbg_ids.reads.plass.faa"
    conda: "plass.yml"
    shell:'''
    plass assemble {input} {output} tmp
    '''

rule cdhit_plass:
    input: "outputs/nbhd_read_plass/{sample}/{gather_genome}.cdbg_ids.reads.plass.faa"
    output: "outputs/nbhd_read_cdhit/{sample}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    conda: "plass.yml"
    shell:'''
    cd-hit -i {input} -o {output} -c 1
    '''

def aggregate_spacegraphcats_gather_matches_plass(wildcards):
    # checkpoint_output produces the output dir from the checkpoint rule.
    checkpoint_output = checkpoints.spacegraphcats_gather_matches.get(**wildcards).output[0]    
    file_names = expand("outputs/nbhd_read_cdhit/{sample}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa",
                        sample = SAMPLES, 
                        gather_genome = glob_wildcards(os.path.join(checkpoint_output, "{gather_genome}.fna.cdbg_ids.reads.fa.gz")).gather_genome)
    return file_names

rule aggregate_spacegraphcats_gather_matches_plass:
    input: aggregate_spacegraphcats_gather_matches_plass
    output: "aggregated_checkpoints/aggregate_spacegraphcats_gather_matches_plass.txt"
    shell:'''
    touch {output}
    '''
rule paladin_index_plass:
    input: "outputs/nbhd_read_cdhit/{sample}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    output: "outputs/nbhd_read_cdhit/{sample}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa.bwt"
    conda: "plass.yml"
    shell: '''
    paladin index -r3 {input}
    '''

rule paladin_align_plass:
    input:
        indx="outputs/nbhd_read_cdhit/{sample}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa.bwt",
        reads="outputs/sgc_genome_queries/{sample}_k31_r1_search_oh0/{gather_genome}.fna.cdbg_ids.reads.fa.gz"
    output: "outputs/nbhd_read_paladin/{sample}/{gather_genome}.sam"
    params: indx = "outputs/nbhd_read_cdhit/{sample}/{gather_genome}.cdbg_ids.reads.plass.cdhit.faa"
    conda: "plass.yml"
    shell:'''
    paladin align -f 125 -t 2 {params.indx} {input.reads} > {output}
    '''


########################################
## PCoA
########################################

rule compare_signatures_cosine:
    input: 
        expand("outputs/filt_sigs_named/{sample}_filt_named.sig", sample = SAMPLES),
    output: "outputs/comp/all_filt_comp.csv"
    conda: "sourmash.yml"
    shell:'''
    sourmash compare -k 31 -p 8 --csv {output} {input}
    '''

rule compare_signatures_jaccard:
    input: 
        expand("outputs/filt_sigs_named/{sample}_filt_named.sig", sample = SAMPLES),
    output: "outputs/comp/all_filt_comp_jaccard.csv"
    conda: "sourmash.yml"
    shell:'''
    sourmash compare --ignore-abundance -k 31 -p 8 --csv {output} {input}
    '''

#rule permanova:

#rule plot_comp:

########################################
## Differential abundance
########################################

rule hash_table_long_unnormalized:
    """
    Unlike the hashtable that is input into the random forest analysis, this
    hash table is not normalized by number of hashes in the filtered signature. 
    Differential expression software that we will be using to calculate 
    differential abundance expects unnormalized counts.
    """
    input: 
        expand("outputs/filt_sigs_named_csv_hmp/{sample}_filt_named.csv", sample = SAMPLES)
    output: csv = "outputs/hash_tables/hmp_unnormalized_abund_hashes_long.csv"
    conda: 'r.yml'
    script: "scripts/all_unnormalized_hash_abund_long.R"
        
rule hash_table_wide_unnormalized:
    input: "outputs/hash_tables/hmp_unnormalized_abund_hashes_long.csv"
    output: "outputs/hash_tables/hmp_unnormalized_abund_hashes_wide.feather"
    run:
        import pandas as pd
        import feather

        ibd = pd.read_csv(str(input), dtype = {"minhash" : "int64", "abund" : "float64", "sample" : "object"})
        ibd_wide=ibd.pivot(index='sample', columns='minhash', values='abund')
        ibd_wide = ibd_wide.fillna(0)
        ibd_wide['sample'] = ibd_wide.index
        ibd_wide = ibd_wide.reset_index(drop=True)
        ibd_wide.columns = ibd_wide.columns.astype(str)
        ibd_wide.to_feather(str(output)) 

#rule differential_abundance_all:
#    input: "outputs/hash_tables/all_unnormalized_abund_hashes_wide.feather"


