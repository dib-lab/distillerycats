import pandas as pd
import feather
from sourmash import signature
import glob
import os
from collections import Counter

metadata_file = "inputs/test_metadata.csv"
m = pd.read_csv(metadata_file, header = 0)
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

#rule download_fastq_files:
#    output: "inputs/raw/{sample}.fastq.gz",
#    run:
#        row = m.loc[m['sample'] == wildcards.sample]
#        fastq = row['download'].values
#        fastq = fastq[0]
#        shell("wget -O {output} {fastq}")

#rule fastp:
#    conda: 'envs/fastp.yml'
#    shell:'''
#    fastp -i in.R1.fq.gz -I in.R2.fq.gz -o out.R1.fq.gz -O out.R2.fq.gz -q 4 -j json_output -l 31 -c
#    '''

rule download_human_db:
    output: "inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz"
    shell:'''
    wget -O {output} https://osf.io/84d59/download
    '''

rule remove_host:
# http://seqanswers.com/forums/archive/index.php/t-42552.html
# https://drive.google.com/file/d/0B3llHR93L14wd0pSSnFULUlhcUk/edit?usp=sharing
    output:
        r1 = 'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        r2 = 'outputs/bbduk/{sample}_R2.nohost.fq.gz',
        human_r1='outputs/bbduk/{sample}_R1.human.fq.gz',
        human_r2='outputs/bbduk/{sample}_R2.human.fq.gz'
    input: 
        r1 = 'outputs/fastp/{sample}_R1.fastp.fq.gz',
        r2 = 'outputs/fastp/{sample}_R2.fastp.fq.gz',
        human='inputs/host/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz'
    conda: 'envs/bbmap.yml'
    shell:'''
    bbduk.sh -Xmx64g t=3 in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} outm={output.human_r1} outm2={output.human_r2} k=31 ref={input.human}
    '''

rule kmer_trim_reads:
    input: 
        'outputs/bbduk/{sample}_R1.nohost.fq.gz',
        'outputs/bbduk/{sample}_R2.nohost.fq.gz'
    output: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    conda: 'envs/sourmash.yml'
    shell:'''
    interleave-reads.py {input} | trim-low-abund.py --gzip -C 3 -Z 18 -M 60e9 -V - -o {output}
    '''

rule compute_signatures:
    input: "outputs/abundtrim/{sample}.abundtrim.fq.gz"
    output: "outputs/sigs/{sample}.sig"
    conda: 'envs/sourmash.yml'
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
    conda: 'envs/sourmash.yml'
    shell:'''
    python scripts/hashvals-to-signature.py -o {output} -k 31 --scaled 2000 --name greater_than_one_count_hashes --filename {input} {input}
    '''

rule filter_signatures_to_greater_than_1_hashes:
    input:
        filt_sig = "outputs/filt_sig_hashes/greater_than_one_count_hashes.sig",
        sigs = "outputs/sigs/{sample}.sig"
    output: "outputs/filt_sigs/{sample}_filt.sig"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig intersect -o {output} -A {input.sigs} -k 31 {input.sigs} {input.filt_sig}
    '''

rule name_filtered_sigs:
    input: "outputs/filt_sigs/{sample}_filt.sig"
    output: "outputs/filt_sigs_named/{sample}_filt_named.sig"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash sig rename -o {output} -k 31 {input} {wildcards.sample}_filt
    '''

rule describe_filtered_sigs:
    input: expand("outputs/filt_sigs_named/{sample}_filt_named.sig", sample = SAMPLES)
    output: "outputs/filt_sigs_named/sig_describe_filt_named_sig.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash signature describe --csv {output} {input}
    '''

rule convert_greater_than_1_signatures_to_csv:
    input: "outputs/filt_sigs_named/{sample}_filt_named.sig"
    output: "outputs/filt_sigs_named_csv/{sample}_filt_named.csv"
    conda: 'envs/sourmash.yml'
    shell:'''
    python scripts/sig_to_csv.py {input} {output}
    '''

rule make_hash_abund_table_long_normalized:
    input: 
        expand("outputs/filt_sigs_named_csv/{sample}_filt_named.csv", sample = SAMPLES)
    output: csv = "outputs/hash_tables/normalized_abund_hashes_long.csv"
    conda: 'envs/r.yml'
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
    conda: 'envs/rf.yml'
    script: "scripts/install_pomona.R"

rule vita_var_sel_rf:
    input:
        info = metadata_file, 
        feather = "outputs/hash_tables/normalized_abund_hashes_wide.feather",
        pomona = "outputs/vita_rf/pomona_install.txt"
    output:
        vita_rf = "outputs/vita_rf/{study}_vita_rf.RDS",
        vita_vars = "outputs/vita_rf/{study}_vita_vars.txt",
        var_filt = "outputs/vita_rf/{study}_ibd_filt.csv"
    params: 
        threads = 32,
        validation_study = "{study}"
    conda: 'envs/rf.yml'
    script: "scripts/vita_rf.R"

rule loo_validation:
    input: 
        var_filt = 'outputs/vita_rf/{study}_var_filt.csv',
        info = metadata_file,
        eval_model = 'scripts/function_evaluate_model.R',
        ggconfusion = 'scripts/ggplotConfusionMatrix.R'
    output: 
        recommended_pars = 'outputs/optimal_rf/{study}_rec_pars.tsv',
        optimal_rf = 'outputs/optimal_rf/{study}_optimal_rf.RDS',
        training_accuracy = 'outputs/optimal_rf/{study}_training_acc.csv',
        training_confusion = 'outputs/optimal_rf/{study}_training_confusion.pdf',
        validation_accuracy = 'outputs/optimal_rf/{study}_validation_acc.csv',
        validation_confusion = 'outputs/optimal_rf/{study}_validation_confusion.pdf'
    params:
        threads = 20,
        validation_study = "{study}"
    conda: 'envs/tuneranger.yml'
    script: "scripts/tune_rf.R"


############################################
## Predictive hash characterization - gather
############################################

rule convert_vita_vars_to_sig:
    input: "outputs/vita_rf/{study}_vita_vars.txt"
    output: "outputs/vita_rf/{study}_vita_vars.sig"
    conda: "envs/sourmash.yml"
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

rule download_gather_refseq:
    output: "inputs/gather_databases/refseq-d2-k31.tar.gz"
    shell:'''
    wget -O {output} https://s3-us-west-2.amazonaws.com/sourmash-databases/2018-03-29/refseq-d2-k31.tar.gz
    '''

rule untar_refseq:
    output: "inputs/gather_databases/refseq-d2-k31.sbt.json"
    input:  "inputs/gather_databases/refseq-d2-k31.tar.gz"
    params: outdir = "inputs/gather_databases"
    shell: '''
    tar xf {input} -C {params.outdir}
    '''

rule gather_vita_vars_all:
    input:
        sig="outputs/vita_rf/{study}_vita_vars.sig",
        db1="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db3="inputs/gather_databases/nayfach-k31.sbt.json",
        db4="inputs/gather_databases/pasolli-mags-k31.sbt.json"
    output: 
        csv="outputs/gather/{study}_vita_vars_all.csv",
        matches="outputs/gather/{study}_vita_vars_all.matches",
        un="outputs/gather/{study}_vita_vars_all.un"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db1} {input.db4} {input.db3} {input.db2}
    '''

rule gather_vita_vars_genbank:
    input:
        sig="outputs/vita_rf/{study}_vita_vars.sig",
        db="inputs/gather_databases/genbank-d2-k31.sbt.json",
    output: 
        csv="outputs/gather/{study}_vita_vars_genbank.csv",
        matches="outputs/gather/{study}_vita_vars_genbank.matches",
        un="outputs/gather/{study}_vita_vars_genbank.un"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db}
    '''

rule gather_vita_vars_refseq:
    input:
        sig="outputs/vita_rf/{study}_vita_vars.sig",
        db="inputs/gather_databases/refseq-d2-k31.sbt.json",
    output: 
        csv="outputs/gather/{study}_vita_vars_refseq.csv",
        matches="outputs/gather/{study}_vita_vars_refseq.matches",
        un="outputs/gather/{study}_vita_vars_refseq.un"
    conda: 'envs/sourmash.yml'
    shell:'''
    sourmash gather -o {output.csv} --save-matches {output.matches} --output-unassigned {output.un} --scaled 2000 -k 31 {input.sig} {input.db}
    '''

rule merge_vita_vars_sig_all:
    input: expand("outputs/vita_rf/{study}_vita_vars.sig", study = STUDY)
    output: "outputs/vita_rf/vita_vars_merged.sig"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash sig merge -o {output} {input}
    '''

rule combine_gather_vita_vars_all:
    output: "outputs/gather/vita_vars_all.csv"
    input: expand("outputs/gather/{study}_vita_vars_all.csv", study = STUDY)
    run:
        import pandas as pd
        
        li = []
        for filename in input:
            df = pd.read_csv(str(filename), index_col=None, header=0)
            df["study"] = str(filename)
            li.append(df)

        frame = pd.concat(li, axis=0, ignore_index=True)
        frame.to_csv(str(output))


checkpoint collect_gather_vita_vars_all_sig_matches:
    input:
        db1="inputs/gather_databases/almeida-mags-k31.sbt.json",
        db2="inputs/gather_databases/genbank-d2-k31.sbt.json",
        db3="inputs/gather_databases/nayfach-k31.sbt.json",
        db4="inputs/gather_databases/pasolli-mags-k31.sbt.json",
        csv="outputs/gather/vita_vars_all.csv"
    output: directory("outputs/gather_matches/")
    run:
        from sourmash import signature
        import pandas as pd

        # load gather results
        df = pd.read_csv(input.csv)

        # for each row, determine which database the result came from
        for index, row in df.iterrows():
            if row["filename"] == "inputs/gather_databases/almeida-mags-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.almeida-mags-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/genbank-d2-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.genbank-d2-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/nayfach-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.nayfach-k31/" + row["md5"]
            elif row["filename"] == "inputs/gather_databases/pasolli-mags-k31.sbt.json":
                sigfp = "inputs/gather_databases/.sbt.pasolli-mags-k31/" + row["md5"]
            # open the signature, parse its name, and write the signature out to a new
            # folder
            sigfp = open(sigfp, 'rt')
            sig = signature.load_one_signature(sigfp)
            out_sig = str(sig.name())
            out_sig = out_sig.split('/')[-1]
            out_sig = out_sig.split(" ")[0]
            out_sig = "outputs/gather_matches/" + out_sig + ".sig"
            with open(str(out_sig), 'wt') as fp:
                signature.save_signatures([sig], fp)      


def aggregate_collect_gather_vita_vars_all_sig_matches(wildcards):
    checkpoint_output = checkpoints.collect_gather_vita_vars_all_sig_matches.get(**wildcards).output[0]  
    file_names = expand("outputs/gather_matches/{genome}.sig", 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}.sig")).genome)
    return file_names
    

rule create_hash_genome_map_gather_vita_vars_all:
    input:
        #genomes = "outputs/gather_matches/{genome}.sig",
        genomes = aggregate_collect_gather_vita_vars_all_sig_matches,
        vita_vars = "outputs/vita_rf/vita_vars_merged.sig"
    output: "outputs/gather_matches_hash_map/hash_to_genome_map_gather_all.csv"
    run:
        from sourmash import signature
        import pandas as pd
        
        sigs = input.genomes
        # read in all genome signatures that had gather 
        # matches for the var imp hashes create a dictionary, 
        # where the key is the genome and the values are the minhashes
        genome_dict = {}
        for sig in sigs:
            sigfp = open(sig, 'rt')
            siglist = list(signature.load_signatures(sigfp))
            loaded_sig = siglist[0] 
            mins = loaded_sig.minhash.get_mins() # Get the minhashes 
            genome_dict[sig] = mins

        # read in vita variables
        sigfp = open(str(input.vita_vars), 'rt')
        vita_vars = sig = signature.load_one_signature(sigfp)
        vita_vars = vita_vars.minhash.get_mins() 

        # generate a list of all minhashes from all genomes
        all_mins = []
        for file in sigs:
            if os.path.getsize(file) > 0:
                sigfp = open(file, 'rt')
                siglist = list(signature.load_signatures(sigfp))
                loaded_sig = siglist[0]
                mins = loaded_sig.minhash.get_mins() # Get the minhashes 
                all_mins += mins

        # define a function where if a hash is a value, 
        # return all key for which it is a value
        def get_all_keys_if_value(dictionary, hash_query):
            genomes = list()
            for genome, v in dictionary.items():
                if hash_query in v:
                    genomes.append(genome)
            return genomes

        # create a dictionary where each vita_vars hash is a key, 
        # and values are the genome signatures in which that hash
        # is contained
        vita_hash_dict = {}
        for hashy in vita_vars:
            keys = get_all_keys_if_value(genome_dict, hashy)
            vita_hash_dict[hashy] = keys

        # transform this dictionary into a dataframe and format the info nicely
        df = pd.DataFrame(list(vita_hash_dict.values()), index = vita_hash_dict.keys())
        df = df.reset_index()
        df = pd.melt(df, id_vars=['index'], var_name= "drop", value_name='genome')
        # remove tmp col drop
        df = df.drop('drop', 1)
        # drop duplicate rows in the df
        df = df.drop_duplicates()
        # write the dataframe to csv
        df.to_csv(str(output), index = False) 


rule download_sourmash_lca_db:
    output: "inputs/gather_databases/gtdb-release89-k31.lca.json.gz"
    shell:'''
    wget -O {output} https://osf.io/gs29b/download
    '''

rule sourmash_lca_classify_vita_vars_all_sig_matches:
    input:
        db = "inputs/gather_databases/gtdb-release89-k31.lca.json.gz",
        genomes = "outputs/gather_matches/{genome}.sig"
    output: "outputs/gather_matches_lca_classify/{genome}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash lca classify --db {input.db} --query {input.genomes} -o {output}
    '''

def aggregate_collect_gather_vita_vars_all_sig_matches_lca_classify(wildcards):
    checkpoint_output = checkpoints.collect_gather_vita_vars_all_sig_matches.get(**wildcards).output[0]  
    file_names = expand("outputs/gather_matches_lca_classify/{genome}.csv", 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}.sig")).genome)
    return file_names

    
rule finished_collect_gather_vita_vars_all_sig_matches_lca_classify:
    input: aggregate_collect_gather_vita_vars_all_sig_matches_lca_classify
    output: "aggregated_checkpoints/finished_collect_gather_vita_vars_all_sig_matches_lca_classify.txt"
    shell:'''
    touch {output}
    '''

rule sourmash_lca_summarize_vita_vars_all_sig_matches:
    input:
        db = "inputs/gather_databases/gtdb-release89-k31.lca.json.gz",
        genomes = "outputs/gather_matches/{genome}.sig"
    output: "outputs/gather_matches_lca_summarize/{genome}.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash lca summarize --db {input.db} --query {input.genomes} -o {output}
    '''

def aggregate_collect_gather_vita_vars_all_sig_matches_lca_summarize(wildcards):
    checkpoint_output = checkpoints.collect_gather_vita_vars_all_sig_matches.get(**wildcards).output[0]  
    file_names = expand("outputs/gather_matches_lca_summarize/{genome}.csv", 
                        genome = glob_wildcards(os.path.join(checkpoint_output, "{genome}.sig")).genome)
    return file_names

    
rule finished_collect_gather_vita_vars_all_sig_matches_lca_summarize:
    input: aggregate_collect_gather_vita_vars_all_sig_matches_lca_summarize
    output: "aggregated_checkpoints/finished_collect_gather_vita_vars_all_sig_matches_lca_summarize.txt"
    shell:'''
    touch {output}
    '''

###################################################
# Predictive hash characterization -- shared hashes
###################################################

########################################
## PCoA
########################################

rule compare_signatures_cosine:
    input: 
        expand("outputs/filt_sigs_named/{sample}_filt_named.sig", sample = SAMPLES),
    output: "outputs/comp/all_filt_comp_cosine.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compare -k 31 -p 8 --csv {output} {input}
    '''

rule compare_signatures_jaccard:
    input: 
        expand("outputs/filt_sigs_named/{sample}_filt_named.sig", sample = SAMPLES),
    output: "outputs/comp/all_filt_comp_jaccard.csv"
    conda: "envs/sourmash.yml"
    shell:'''
    sourmash compare --ignore-abundance -k 31 -p 8 --csv {output} {input}
    '''

rule permanova_jaccard:
    input: 
        comp = "outputs/comp/all_filt_comp_jaccard.csv",
        info = metadata_file,
        sig_info = "outputs/filt_sigs_named/sig_describe_filt_named_sig.csv"
    output: 
        perm = "outputs/comp/all_filt_permanova_jaccard.csv"
    conda: "envs/vegan.yml"
    script: "scripts/run_permanova.R"

rule permanova_cosine:
    input: 
        comp = "outputs/comp/all_filt_comp_cosine.csv",
        info = metadata_file,
        sig_info = "outputs/filt_sigs_named/sig_describe_filt_named_sig.csv"
    output: 
        perm = "outputs/comp/all_filt_permanova_cosine.csv"
    conda: "envs/vegan.yml"
    script: "scripts/run_permanova.R"

rule plot_comp_jaccard:
    input:
        comp = "outputs/comp/all_filt_comp_jaccard.csv",
        info = metadata_file
    output: 
        study = "outputs/comp/study_plt_all_filt_jaccard.pdf",
        var = "outputs/comp/var_plt_all_filt_jaccard.pdf"
    conda: "envs/ggplot.yml"
    script: "scripts/plot_comp.R"

rule plot_comp_cosine:
    input:
        comp = "outputs/comp/all_filt_comp_cosine.csv",
        info = metadata_file
    output: 
        study = "outputs/comp/study_plt_all_filt_cosine.pdf",
        var = "outputs/comp/var_plt_all_filt_cosine.pdf"
    conda: "envs/ggplot.yml"
    script: "scripts/plot_comp.R"

