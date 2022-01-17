#!/usr/bin/env python3

# convert large parquet files into directories with smaller partition sizes and snappy compression

from pathlib import Path
import sys
#from ratschlab_common import spark_wrapper
import spark_wrapper

def main():
    if len(sys.argv) != 5:
        sys.exit("need 4 arguments")

    input_dir = Path(sys.argv[1])

    output_dir = Path(sys.argv[2])

    cores = int(sys.argv[3])
    mem = int(sys.argv[4])

    output_dir.mkdir(exist_ok=True, parents=True)

    spark_cfg = spark_wrapper.default_spark_config(cores, mem,
                                                   extra_java_options='-XX:+UseG1GC')

    with spark_wrapper.create_spark_session_from_config(spark_cfg) as spark:
        for f in input_dir.glob('*.pq.*'):
            print(f"Working on {str(f)}")
            output = output_dir/f.name

#            if f.name != 'SRR665538.aligned_germline_graph_kmer_SegmExpr.pq.gz':
#                # skipping corrupt file
#                print('skipping {}'.format(f.name))
#                continue
#
            if (output/'_SUCCESS').exists():
                print(f"output {output} already exists")
                continue

            df_in = spark.read.parquet(str(f))
            df_in.write.parquet(str(output_dir/f.name), mode='overwrite')


if __name__ == "__main__":
    main()
