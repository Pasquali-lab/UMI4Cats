#!/usr/bin/env python
# coding: utf-8

# In[ ]:


####################################################################################################################

# UMI4Cats python functions
# description in R

####################################################################################################################

def prep(raw_dir,
        wk_dir,
        bait_seq,
        bait_pad,
        res_e,
        fastqmultx):

    # import modules
    import re
    import os
    from itertools import islice
    import glob
    import datetime
    import subprocess
    from art import tprint
    import sys

    # create directory
    prepDir = os.path.join(wk_dir, 'prep')

    try:
        os.mkdir(prepDir)
    except:
        pass

    # Starting message
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

    def catPrint():
        print('                     _________')
        print('                    /         \\')
        print('                  | ___________ |')
        print('      /\\_/\\       | |         | |')
        print(' /\\  / o o \\      | |         | |')
        print('//\\\\ \\~(*)~/      | |_________| |')
        print('`  \\/   ^ /       \\_____________/ ')
        print('   | \\|| ||       / \"\"\"\"\"\"\"\"\"\"\" \\')
        print('   \\ \'|| ||      / ::::::::::::: \\')
        print('    \\)()-())    (_________________)')

    tprint("UMI4Cats prep")

    catPrint()

    print('\n')

    print(bcolors.UNDERLINE + 'Preprocessing' + bcolors.ENDC + '\n\n',
            'Raw directory:', raw_dir, '\n',
            'Work directory:', wk_dir, '\n',
            'Bait sequence:', bait_seq, '\n',
            'Bait pad:', bait_pad, '\n',
            'Restriction enzyme sequence:', res_e, '\n')

    print('\n')

    # necessary for showing output using reticulate in R
    sys.stdout.flush()

    # descompress fastq.gz if it is necessary
    compressedFiles = glob.glob(os.path.join(raw_dir, '*.gz'))

    if len(compressedFiles) != 0:
        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Descompressing', compressedFiles)
        sys.stdout.flush()
        ['gunzip -c {file} > {decompFile}'.format(file = file, decompFile = os.path.splitext(file)[0]) for file in compressedFiles]

        for file in compressedFiles:
            cmd = 'gunzip -k -c {file} > {decompFile}'.format(file = file, decompFile = os.path.splitext(file)[0])
            subprocess.call(cmd, shell='True')
            print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
                  'Finished descompression of', compressedFiles)
            sys.stdout.flush()
    else:
        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Non files to decompress')
        sys.stdout.flush()
    # insert umi identifier (10 first bp of R2) to header of both R1 R2 files


    R2Files = glob.glob(os.path.join(raw_dir, '*_R2.fq')) + glob.glob(os.path.join(raw_dir, '*_R2.fastq'))
    R1Files = glob.glob(os.path.join(raw_dir, '*_R1.fq')) + glob.glob(os.path.join(raw_dir, '*_R1.fastq'))

    for R1, R2 in zip(R1Files, R2Files):
        nameFile = re.sub("_R1.fastq", "", os.path.basename(R1))

        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Inserting UMI identifier in', nameFile)
        sys.stdout.flush()
        umiFile1 = os.path.join(prepDir, nameFile + '_umi_R1.fastq')
        umiFile2 = os.path.join(prepDir, nameFile + '_umi_R2.fastq')

        cmd = 'cat {R2} | paste - - - - | '                 'awk \'BEGIN {{FS="\\t"; OFS="\\t"}} {{split($1, a, " ")}}; {{print a[1]":"substr($2, 1, 10), $2, $3, $4}}\' |'                 'sed \'s/\\t/\\n/g\' > {umiFile2}'.format(R2 = R2, umiFile2 = umiFile2)

        subprocess.call(cmd, shell='True')

        cmd = 'paste '                 '<(cat {umiFile2} | paste - - - - | awk \'BEGIN {{FS="\\t"}} {{print $1}}\') '                 '<(cat {R1} | paste - - - - | awk \'BEGIN {{FS="\\t"; OFS="\\t"}} {{print $2, $3, $4}}\') | '                 'sed \'s/\\t/\\n/g\' > {umiFile1}'.format(umiFile2 = umiFile2, R1 = R1, umiFile1 = umiFile1)

        subprocess.call(cmd, shell='True', executable='/bin/bash')

        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Finished insertion UMI identifier in', nameFile)
        sys.stdout.flush()
    # Filter reads that not present bait in R1

    # prepare barcode file
    barcodes = os.path.join(prepDir, nameFile + "_barcodes.txt")

    # in some cases one of the bait elements it is not present
    bait = [bait_seq, bait_pad, res_e]

    with open (barcodes, "w") as barcodeFile:
        barcodeFile.write("prefiltered\t" + "".join(bait))

    # filter using fastq-multx
    umisR1 = glob.glob(os.path.join(prepDir, '*umi_R1.fastq'))
    umisR2 = glob.glob(os.path.join(prepDir, '*umi_R2.fastq'))

    for umiR1, umiR2 in zip(umisR1, umisR2):

        nameFile = re.sub("_umi_R1.fastq", "", os.path.basename(umiR1))

        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Filtering', nameFile)
        sys.stdout.flush()

        prefiltered1 = os.path.join(prepDir, nameFile + '_%_R1.fastq')
        prefiltered2 = os.path.join(prepDir, nameFile + '_%_R2.fastq')

        cmd = '{fastqmultx} '                 '-x '                 '-m 0 '                 '-b {barcodes} '                 '{umiR1} {umiR2} '                 '-o {prefiltered1} '                 '-o {prefiltered2} '.format(fastqmultx = fastqmultx,
                                                barcodes = barcodes,
                                                umiR1 = umiR1,
                                                umiR2 = umiR2,
                                                prefiltered1 = prefiltered1,
                                                prefiltered2 = prefiltered2)

        subprocess.call(cmd, shell='True')

        # remove annoying header modification done by fastq-multx and r1 or r2 is added for each file respectively

        prefiltered1 = os.path.join(prepDir, nameFile + '_prefiltered_R1.fastq')
        prefiltered2 = os.path.join(prepDir, nameFile + '_prefiltered_R2.fastq')
        filtered1 = os.path.join(prepDir, nameFile + '_filtered_R1.fastq')
        filtered2 = os.path.join(prepDir, nameFile + '_filtered_R2.fastq')

        cmd = 'cat {prefiltered1} | paste - - - - | '                 'awk \'BEGIN {{FS="\\t"; OFS="\\t"}}; {{split($0, a, " ")}}; {{print a[1]":R1", $2, $3, $4}}\' | '                 'sed \'s/\\t/\\n/g\' > {filtered1}'.format(prefiltered1 = prefiltered1,
                                                           filtered1 = filtered1)

        subprocess.call(cmd, shell='True', executable='/bin/bash')

        cmd = 'cat {prefiltered2} | paste - - - - | '                 'awk \'BEGIN {{FS="\\t"; OFS="\\t"}}; {{split($0, a, " ")}}; {{print a[1]":R2", $2, $3, $4}}\' | '                 'sed \'s/\\t/\\n/g\' > {filtered2}'.format(prefiltered2 = prefiltered2,
                                                           filtered2 = filtered2)
        subprocess.call(cmd, shell='True', executable='/bin/bash')

        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Finished filtering of', nameFile)
        sys.stdout.flush()

def split(wk_dir,
          res_e,
          cut_pos):

    # import modules
    import re
    import os
    from itertools import islice
    import glob
    import datetime
    from art import tprint
    import sys

    # create directory
    prepDir = os.path.join(wk_dir, 'prep')
    splitDir = os.path.join(wk_dir, 'split')

    try:
        os.mkdir(splitDir)
    except:
        pass

    # Starting message
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

    def catPrint():
        print('                     _________')
        print('                    /         \\')
        print('                  | ___________ |')
        print('      /\\_/\\       | |         | |')
        print(' /\\  / o o \\      | |         | |')
        print('//\\\\ \\~(*)~/      | |_________| |')
        print('`  \\/   ^ /       \\_____________/ ')
        print('   | \\|| ||       / \"\"\"\"\"\"\"\"\"\"\" \\')
        print('   \\ \'|| ||      / ::::::::::::: \\')
        print('    \\)()-())    (_________________)')

    tprint("UMI4Cats split")

    catPrint()

    print('\n')

    print(bcolors.UNDERLINE + 'Splitting' + bcolors.ENDC + '\n\n',
            'Work directory:', wk_dir, '\n',
            'Restriction enzyme sequence:', res_e, '\n',
            'Restriction enzyme cut position:', cut_pos, '\n',)

    print('\n')
    sys.stdout.flush()

    # define split function
    def splitFastq(input_file, output_file, res_e, cut_pos, Rtype):

        # correct cutting
        if Rtype == 'R2':
            cut_pos = -cut_pos

        with open(input_file, "r") as f_in, open(output_file, "a") as f_out:
            while True:
                fastqLines = list(islice(f_in, 4))
                if not fastqLines:
                    break

                # process fastqLines
                header = fastqLines[0].rstrip()
                sequence = fastqLines[1].rstrip()
                quality = fastqLines[3].rstrip()

                # define where to cut
                splitPoints = [m.start() for m in re.finditer(res_e, sequence)]
                splitPoints = [0] + splitPoints
                splitPoints.append(len(sequence))

                # correct cut by cutting position of the restriction enzyme
                splitPoints = [ x + cut_pos for x in splitPoints ]

                # split by restiction enzyme
                sequences = [sequence[i: j] for i, j in zip(splitPoints, splitPoints[1:])]
                qualities = [quality[i: j] for i, j in zip(splitPoints, splitPoints[1:])]

                for sequence, quality in zip(sequences, qualities):
                    f_out.write(header + "\n" + sequence + "\n" + "+\n" + quality + "\n")

        print('Finished {file}'.format(file = input_file))
        sys.stdout.flush()

    filesTosplitR1 = glob.glob(os.path.join(prepDir, '*_filtered_R1.fastq'))
    filesTosplitR2 = glob.glob(os.path.join(prepDir, '*_filtered_R2.fastq'))

    cut_pos = int(cut_pos)

    # split files
    for input_file_R1, input_file_R2 in zip(filesTosplitR1, filesTosplitR2):
        nameFile = re.sub("_filtered_R1.fastq", "", os.path.basename(input_file_R1))

        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Splitting', nameFile)
        sys.stdout.flush()

        output_file = os.path.join(splitDir, nameFile + '_split.fastq')
        open(output_file, 'w').close()

        splitFastq(input_file_R1, output_file, res_e, cut_pos, 'R1')
        splitFastq(input_file_R2, output_file, res_e, cut_pos, 'R2')

        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Finished splitting of', nameFile)
        sys.stdout.flush()

def alignment(wk_dir,
                threads,
                bowtie2,
                ref_gen,
                samtools,
                bait_seq,
                bait_pad,
                res_e):

    # import modules
    import re
    import os
    import glob
    import datetime
    import subprocess
    from art import tprint
    import sys

    # create directory
    alignmentDir = os.path.join(wk_dir, 'alignment')
    splitDir = os.path.join(wk_dir, 'split')

    try:
        os.mkdir(alignmentDir)
    except:
        pass

    # Starting message
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

    def catPrint():
        print('                     _________')
        print('                    /         \\')
        print('                  | ___________ |')
        print('      /\\_/\\       | |         | |')
        print(' /\\  / o o \\      | |         | |')
        print('//\\\\ \\~(*)~/      | |_________| |')
        print('`  \\/   ^ /       \\_____________/ ')
        print('   | \\|| ||       / \"\"\"\"\"\"\"\"\"\"\" \\')
        print('   \\ \'|| ||      / ::::::::::::: \\')
        print('    \\)()-())    (_________________)')

    tprint("UMI4Cats alignment")

    catPrint()

    print('\n')

    print(bcolors.UNDERLINE + 'Alignment' + bcolors.ENDC + '\n\n',
            'Work directory:', wk_dir, '\n',
            'Bait sequence:', bait_seq, '\n',
            'Bait pad:', bait_pad, '\n',
            'Restriction enzyme sequence: ', res_e, '\n',
            'bowtiw2 path:', bowtie2, '\n'
            'samtools path:', samtools, '\n'
            'Reference genome used:', ref_gen, '\n'
            'Number of threads used:', threads)

    print('\n')
    sys.stdout.flush()

    splitFiles = glob.glob(os.path.join(splitDir, '*_split.fastq'))

    # get coordinates of viewpoint using bowtie2
    viewpoint = bait_seq + bait_pad + res_e
    index = os.path.splitext(ref_gen)[0]

    cmd = "bowtie2 --quiet -x {index} -c {viewpoint} -N 0 | "\
          "samtools view | "\
          "awk \'{{print $3,$4}}\'".format(index = index,
                                            viewpoint = viewpoint)

    viewpointPos = subprocess.check_output(cmd, shell=True)
    viewpointPos = viewpointPos.decode("utf-8")
    viewpointPos = re.sub('\n', '', viewpointPos)
    chrVP = viewpointPos.split(' ')[0]
    startVP = viewpointPos.split(' ')[1]
    endVP = int(startVP) + len(viewpoint)

    # align splited files
    for file in splitFiles:
        nameFile = re.sub("_split.fastq", "", os.path.basename(file))
        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Aligning', nameFile)
        sys.stdout.flush()

        sam = os.path.join(alignmentDir, nameFile + '.sam')
        bam = os.path.join(alignmentDir, nameFile + '.bam')
        samFiltered = os.path.join(alignmentDir, nameFile + '_filtered.sam')
        log = os.path.join(alignmentDir, nameFile + '.log')
        ref_gen = os.path.splitext(ref_gen)[0]

        cmd = '{bowtie2} '\
                '-p {threads} '\
                '-x {ref_gen} '\
                '-U {file} '\
                '-S {sam} '\
                '2> {log} '.format(bowtie2 = bowtie2,
                                              threads = threads,
                                              ref_gen = ref_gen,
                                              file = file,
                                              sam = sam,
                                              log = log)

        subprocess.call(cmd, shell=True)

        cmd = 'samtools view {sam} -@ {threads} -Sb  | '\
              'samtools sort -@ {threads} - -o {bam} ; '\
              'samtools index -@ {threads} {bam}'.format(threads = threads,
                                                                  sam = sam,
                                                                  bam = bam)
        subprocess.call(cmd, shell=True)

        # without header for being able to load in pandas
        # only chrom where the view point is present and mapq of 42 at least
        cmd = 'samtools view {bam} -@ {threads} -q 42 {chrVP} > {samFiltered}'.format(threads = threads,
                                                                                    bam = bam,
                                                                                    chrVP = chrVP,
                                                                                    samFiltered = samFiltered)
        subprocess.call(cmd, shell=True)
        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Finished alignment of', nameFile)
        sys.stdout.flush()

def umiCounter(wk_dir,
            bowtie2,
            ref_gen,
            samtools,
            genomic_track,
            bait_seq,
            bait_pad,
            res_e):

    # import modules
    import re
    import os
    import glob
    import datetime
    import subprocess
    import pandas as pd
    import pyranges as pr
    import warnings
    from art import tprint
    import sys

    warnings.filterwarnings('ignore')

    # create directory
    alignmentDir = os.path.join(wk_dir, 'alignment')
    fragsDir = os.path.join(wk_dir, 'frags')

    try:
        os.mkdir(fragsDir)
    except:
        pass

            # Starting message
    class bcolors:
        HEADER = '\033[95m'
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        WARNING = '\033[93m'
        FAIL = '\033[91m'
        ENDC = '\033[0m'
        BOLD = '\033[1m'
        UNDERLINE = '\033[4m'

    def catPrint():
        print('                     _________')
        print('                    /         \\')
        print('                  | ___________ |')
        print('      /\\_/\\       | |         | |')
        print(' /\\  / o o \\      | |         | |')
        print('//\\\\ \\~(*)~/      | |_________| |')
        print('`  \\/   ^ /       \\_____________/ ')
        print('   | \\|| ||       / \"\"\"\"\"\"\"\"\"\"\" \\')
        print('   \\ \'|| ||      / ::::::::::::: \\')
        print('    \\)()-())    (_________________)')

    tprint("UMI4Cats umiCounter")

    catPrint()

    print('\n')

    print(bcolors.UNDERLINE + 'Counting UMIs' + bcolors.ENDC + '\n\n',
            'Work directory:', wk_dir, '\n',
            'Bait sequence:', bait_seq, '\n',
            'Bait pad:', bait_pad, '\n',
            'Restriction enzyme sequence:', res_e, '\n',
            'bowtiw2 path:', bowtie2, '\n'
            'samtools path:', samtools, '\n'
            'Reference genome used:', ref_gen, '\n'
            'Genomic track used:', genomic_track)

    print('\n')
    sys.stdout.flush()

    # get coordinates of viewpoint using bowtie2
    viewpoint = bait_seq + bait_pad + res_e
    index = os.path.splitext(ref_gen)[0]
    cmd = "bowtie2 --quiet -x {index} -c {viewpoint} -N 0 | "\
          "samtools view | "\
          "awk \'{{print $3,$4}}\'".format(index = index,
                                              viewpoint = viewpoint)

    viewpointPos = subprocess.check_output(cmd, shell=True)
    viewpointPos = viewpointPos.decode("utf-8")
    viewpointPos = re.sub('\n', '', viewpointPos)
    chrVP = viewpointPos.split(' ')[0]
    startVP = viewpointPos.split(' ')[1]
    endVP = int(startVP) + len(viewpoint)

    # count UMIs for every ligation

    samFiles = glob.glob(os.path.join(alignmentDir, '*_filtered.sam'))

    for sam in samFiles:

        nameFile = re.sub("_filtered.sam", "", os.path.basename(sam))

        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Genereting UMI counts of', nameFile)
        sys.stdout.flush()

        # create df from sam
        dfSam = pd.read_csv(sam,
                            sep = '\t',
                            header = None,
                            usecols = list(range(5)) + [9])

        samColumns = ['header',
                     'string',
                     'chr',
                     'start',
                     'mapq',
                     'seq'
                     ]
        dfSam.columns = samColumns

        # select only mapped reads and rename forward and reverse reads
        dfSam = dfSam[(dfSam.string == 16) | (dfSam.string == 0)]
        dfSam.string[dfSam.string == 16] = "-"
        dfSam.string[dfSam.string == 0] = "+"

        # define end coor and header
        dfSam['end'] = dfSam.start + dfSam.seq.str.len() -1
        dfSam = dfSam[['header', 'string', 'chr', 'start', 'end', 'mapq', 'seq']]

        # load genomic track and transform to df
        dfGenomicTrack = pd.read_csv(genomic_track, sep = '\t', header = None, index_col = None)
        genomicTrackColumns = ['chr', 'start', 'end']
        dfGenomicTrack.columns = genomicTrackColumns
        dfGenomicTrackChr = dfGenomicTrack[dfGenomicTrack['chr'] == chrVP]

        # correct df sam and create pyranges object
        dfSam = dfSam[['chr', 'start', 'end', 'string', 'header', 'mapq', 'seq']]
        dfSam.columns = ['Chromosome', 'Start', 'End', 'String', 'header', 'mapq', 'seq']
        prSam = pr.PyRanges(dfSam)

        # correct df genomic ranges and create pyranges object
        dfGenomicTrackChr['string'] = '+'
        dfGenomicTrackChr.columns = ['Chromosome', 'Start', 'End', 'String']
        prGenomicTrackChr = pr.PyRanges(dfGenomicTrackChr)

        # overlap segments from sam and frags from genomic ranges
        dfFragSeg = prSam.nearest(prGenomicTrackChr).df

        dfFragSeg.columns = ['Chromosome', 'StartSeg', 'EndSeg',
                                'String', 'header', 'mapq', 'seq',
                                'StartFrag', 'EndFrag', 'StringFrag',
                                'Distance']

        dfFragSeg = dfFragSeg[['Chromosome', 'StartSeg', 'EndSeg',
                                'StartFrag', 'EndFrag', 'String',
                                'header', 'mapq', 'seq']]

        # correct dataframe of overlapping
        dfFragSeg.drop_duplicates(inplace=True)
        dfFragSeg = dfFragSeg.reset_index(drop=True)
        dfFragSeg['end'] = dfFragSeg['header'].str.split(":").apply(lambda x: x[-1])
        newHeader = dfFragSeg['header'].str.split(":").apply(lambda x: x[:-1])
        dfFragSeg['header'] = list(map(lambda x: ":".join(x), newHeader))

        # get viewpoint fragment
        dfViewPoint = pd.DataFrame([chrVP, startVP, endVP]).transpose()
        dfViewPoint.columns = ['Chromosome', 'Start', 'End']
        prViewPoint = pr.PyRanges(dfViewPoint)
        dfFragViewPoint = prViewPoint.nearest(prGenomicTrackChr).df
        dfFragViewPoint = dfFragViewPoint.iloc[:,[0,3,4]]
        dfFragViewPoint.columns = ['Chromosome', 'StartFrag', 'EndFrag']

        # select viewpoint reads and contact reads
        keys = list(dfFragViewPoint.columns.values)
        i1 = dfFragSeg.set_index(keys).index
        i2 = dfFragViewPoint.set_index(keys).index
        dfFragSegViewPoint = dfFragSeg[i1.isin(i2)]
        dfFragSegContacts = dfFragSeg[~i1.isin(i2)]
        dfFragSegViewPoint.reset_index(drop=True, inplace = True)
        dfFragSegContacts.reset_index(drop=True, inplace = True)

        # generate ligations merging contacts with viewpoint using as a reference header of fastq
        dfViewPointContact = pd.merge(dfFragSegViewPoint, dfFragSegContacts, on='header', how='inner')
        newColumns = list(map(lambda x: re.sub("x", "Frag1", x), dfViewPointContact.columns))
        newColumns = list(map(lambda x: re.sub("y", "Frag2", x), newColumns))
        dfViewPointContact.columns = newColumns

        # generate umi and pos columns
        dfViewPointContact['umi'] = dfViewPointContact['header'].apply(lambda x: x.split(':')[-1])
        dfViewPointContact['pos'] = dfViewPointContact['Chromosome_Frag1'].astype(str) +\
                                      "_" +\
                                      dfViewPointContact['StartSeg_Frag1'].astype(str) +\
                                      "_" +\
                                      dfViewPointContact['EndSeg_Frag1'].astype(str) +\
                                      "_" +\
                                      dfViewPointContact['Chromosome_Frag2'].astype(str) +\
                                      "_" +\
                                      dfViewPointContact['StartSeg_Frag2'].astype(str) +\
                                      "_" +\
                                      dfViewPointContact['EndSeg_Frag2'].astype(str)

        newColumns = list(map(lambda x: re.sub("x", "Frag1", x), dfViewPointContact.columns))
        newColumns = list(map(lambda x: re.sub("y", "Frag2", x), newColumns))
        dfViewPointContact.columns = newColumns

        # umi filtering colapsing ligations with the same position or with an umi with a less than 2 mismatches

        # colapse ligations with same pos
        dfColapsedPos = dfViewPointContact.drop_duplicates(subset = 'pos')
        dfColapsedPos.reset_index(drop = True, inplace = True)

        dfColapsedPosSub = dfColapsedPos[['umi', 'pos']]

        # colapse ligations with umis with less than 3 mismatches
        def mismatches(a, b):
            count = sum(c1!=c2 for c1,c2 in zip(a, b))
            return(count)

        dfInput = dfColapsedPosSub
        dfColapsedUmi = pd.DataFrame()

        for index in dfInput.index:
            dfComparison = dfInput.loc[dfInput.index != index, :]
            umi = dfInput['umi'][index]
            dfFiltered = dfComparison[[mismatches(x, umi) < 2 for x in dfComparison['umi']]]
            dfColapsedUmi = pd.concat([dfColapsedUmi, dfFiltered])
            dfColapsedUmi.drop_duplicates(inplace = True)

        dfTrueContacts = dfColapsedPosSub[~dfColapsedPosSub.isin(dfColapsedUmi)].dropna()
        dfTrueContacts = dfTrueContacts.merge(dfColapsedPos, how = 'inner')
        dfTrueContactsClean = dfTrueContacts[['Chromosome_Frag1', 'StartFrag_Frag1',
                                              'Chromosome_Frag2', 'StartFrag_Frag2']]

        # count umis for ligation
        dfCounts = pd.DataFrame(dfTrueContactsClean['StartFrag_Frag2'].value_counts())
        dfCounts['counts'] = dfCounts['StartFrag_Frag2']
        dfCounts['StartFrag_Frag2'] = dfCounts.index
        dfCounts.reset_index(drop = True, inplace = True)
        dfFinalCounts = dfCounts.merge(dfTrueContactsClean, how = 'left')

        # define contact position of viewpoint frag being the end fragment most proximal to the construction
        startDiff = abs(dfFragViewPoint['StartFrag'][0] - int(startVP))

        endDiff = abs(dfFragViewPoint['EndFrag'][0] - int(startVP))

        if startDiff < endDiff:
          StartFrag = dfFragViewPoint['StartFrag'][0]
        elif startDiff > endDiff:
          StartFrag = dfFragViewPoint['EndFrag'][0]

        dfFinalCounts['StartFrag_Frag1'] = StartFrag

        dfFinalCounts = dfFinalCounts[['Chromosome_Frag1', 'StartFrag_Frag1', 'Chromosome_Frag2', 'StartFrag_Frag2', 'counts']]
        dfFinalCounts.drop_duplicates(inplace = True)
        dfFinalCounts.reset_index(drop = True, inplace = True)

        # create directory
        rstDir = os.path.join(wk_dir, 'rst')
        try:
            os.mkdir(rstDir)
        except:
            pass

        # save counts
        outFile = os.path.join(rstDir, nameFile + '_umi_counts.tsv')
        dfFinalCounts.to_csv(outFile, sep = '\t', index = False, header = False)

        print('[' + str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")) + ']',
              'Finished generation of UMI counts of', nameFile, '\n')
        print('Output saved at', outFile)
        sys.stdout.flush()

