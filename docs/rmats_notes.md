## Installation
```
module load samtools 
conda create -n rmats "python==3.6" cython numpy
conda activate rmats

git clone https://github.com/Xinglab/rmats-turbo.git
cd rmats-turbo

# change lines 11, 12 from rMATS_pipeline/setup.py
# os.environ['CC'] = '/usr/bin/cc'
# os.environ['CXX'] = '/usr/bin/c++'

./build_rmats --no-paired-model > build.log

# comment lines, 39, 40 and 52 from test_rmats
./test_rmats
# three test will fail
# test (tests.paired_stats.test.FilteredEventTest) ... FAIL
# test (tests.paired_stats.test.OneEventTest) ... FAIL
# test (tests.paired_stats.test.TwoEventTest) ... FAIL
```
## Notes on the method

- For multi gene locus, rMATS may mismatch the origin of the gene with novel AS event
https://github.com/Xinglab/rmats-turbo/issues/79#issuecomment-763148314

- rMATS only check annotated introns for IR
https://github.com/Xinglab/rmats-turbo/issues/65#issuecomment-740722530
