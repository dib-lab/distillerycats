""" link or download databases """

import os, sys

database_dir = config['database_dir']
db_logs = os.path.join(database_dir, "logs")

database_info = config["available_databases"]

# if a url, download. Otherwise, link.
urls_begin = ["http", "ftp"]

localrules: get_sbt

rule get_sbt:
    output:
        os.path.join(database_dir, "{database}.sbt.zip")
    params:
        sbt_info= lambda w: database_info.at[w.database,"path"]
    log: os.path.join(db_logs, "get_dbs", "{database}.sbt.get")
    threads: 1
    resources:
        mem_mb=1000,
        runtime=600
    run:
        if params.sbt_info.startswith(tuple(urls_begin)):
            shell("curl -L {params.sbt_info}  > {output}")
        else:
            full_input = os.path.abspath(str(params.sbt_info))
            full_output = os.path.abspath(str(output))
            shell("ln -s {full_input} {full_output} 2> {log}")

localrules: get_dbinfo

rule get_dbinfo:
    output:
        os.path.join(database_dir, "{database}.info.csv")
    params:
         csv_info= lambda w: database_info.at[w.database, 'info_path']
    log: os.path.join(db_logs, "get_dbs", "{db_basename}.info.get")
    threads: 1
    resources:
        mem_mb=1000,
        runtime=600
    run:
        if params.csv_info.startswith(tuple(urls_begin)):
            shell("curl -L {params.csv_info}  > {output}")
        else:
            full_input = os.path.abspath(str(params.csv_info))
            full_output = os.path.abspath(str(output))
            shell("ln -s {full_input} {full_output} 2> {log}")

