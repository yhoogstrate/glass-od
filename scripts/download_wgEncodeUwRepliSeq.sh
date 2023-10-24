#!/usr/bin/env sh

wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig
wget -P data/wgEncodeUwRepliSeq/  http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig





CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBg02esWaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqBjWaveSignalRep2.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm06990WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12801WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12812WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12813WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqGm12878WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHelas3WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHepg2WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqHuvecWaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqImr90WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqK562WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqMcf7WaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqNhekWaveSignalRep1.bigWig.hg38.out
CrossMap.py bigwig data/reference/hg19ToHg38.over.chain.gz data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig data/wgEncodeUwRepliSeq/wgEncodeUwRepliSeqSknshWaveSignalRep1.bigWig.hg38.out



