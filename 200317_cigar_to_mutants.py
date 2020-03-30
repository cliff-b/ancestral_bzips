import re

gencode = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W'}


def revcomp(dna):
    trans = str.maketrans('ACTG', 'TGAC')
    complement = dna.translate(trans)
    return complement[::-1]


def samreader(filename):
    cigardic = {}
    with open(filename) as f:
        for line in f:
            linelist = line.split(sep="\t")
            if linelist[0] not in cigardic:
                cigardic[linelist[0]] = {'name': linelist[2][:-5], 'cigar': linelist[5], 'sequence': linelist[9]}
    return cigardic


def fastareader(filename, side):
    fastadic = {}
    with open(filename) as f:
        for line in f:
            line = line.strip()
            if line[-3] == '1' and line[-1] == side:
                seq = f.readline().strip()
                fastadic[line[1:-4]] = seq
    return fastadic


def cigarparse(cigardic):
    parsedic = {}
    for bc in cigardic:
        curcig = cigardic[bc]['cigar']
        matches = re.findall(r'(\d+)([A-Z=]{1})', curcig)
        parsecig = [(int(m[0]), m[1]) for m in matches]
        startpoint = 0
        startlist = []
        for point in parsecig:
            if startpoint < 126:
                startpoint += point[0]
                startlist.append(point)
        startlist[-1] = (startlist[-1][0] - (startpoint - 126), startlist[-1][1])  # this gets us to the start of the reading frame for the x peptide
        if all(a[1] == '=' or a[1] == 'X' for a in startlist):
            if any(a[1] == 'X' for a in startlist):
                parsedic[bc] = {'name': cigardic[bc]['name'], 'cigar': cigardic[bc]['cigar'], 'sequence': cigardic[bc]['sequence']}
            fromstart = 0
            for i in reversed(startlist):
                if i[1] == 'X':
                    startframe = fromstart % 3
                    endframe = (126 - fromstart - i[0]) % 3
                    inframerevcodons = cigardic[bc]['sequence'][126-fromstart-i[0]-endframe:126-fromstart+startframe]
                    # print("barcode ", bc, " with cigar string ", startlist, " and codons ", inframerevcodons, " fromstart ", fromstart, " startframe ", startframe, " endframe ", endframe )
                    for j in range(0, len(inframerevcodons), 3):
                        try:
                            parsedic[bc]['codons'].append(inframerevcodons[j:j+3])
                        except KeyError:
                            parsedic[bc]['codons'] = [inframerevcodons[j:j+3]]
                        try:
                            parsedic[bc]['positions'].append((126-fromstart-i[0]-endframe + j, 126-fromstart-i[0]-endframe + j + 3))
                        except KeyError:
                            parsedic[bc]['positions'] = [(126-fromstart-i[0]-endframe + j, 126-fromstart-i[0]-endframe + j + 3)]
                fromstart += i[0]
    return parsedic


def cigarparse_y(cigardic):
    parsedic = {}
    for bc in cigardic:
        if len(cigardic[bc]['sequence']) > 164:
            curcig = cigardic[bc]['cigar']
            matches = re.findall(r'(\d+)([A-Z=]{1})', curcig)
            parsecig = [(int(m[0]), m[1]) for m in matches]
            startpoint, endpoint = 0, 0
            startlist = []
            done = False
            for point in reversed(parsecig):
                if startpoint < 164:
                    if endpoint + point[0] < 35 and not done:
                        endpoint += point[0]
                        startpoint += point[0]
                    else:
                        done = True
                        startpoint += point[0]
                        startlist.append(point)
            startlist[-1] = (startlist[-1][0] - (startpoint - 164), startlist[-1][1])  # this gets us to the start of the reading frame for the y peptide
            startlist[0] =(startlist[0][0] - (35 - endpoint), startlist[0][1])
            if startlist[0][0] == 0:
                startlist = startlist[1:]
            if all(a[1] == '=' or a[1] == 'X' for a in startlist):
                if any(a[1] == 'X' for a in startlist):
                    parsedic[bc] = {'name': cigardic[bc]['name'], 'cigar': cigardic[bc]['cigar'], 'sequence': cigardic[bc]['sequence']}
                fromend = 0
                for i in startlist:
                    if i[1] == 'X':
                        inframerevcodons = cigardic[bc]['sequence'][-(35 + (fromend + i[0]+2)//3*3):-(35 + (fromend//3)*3)]
                        for j in range(0, len(inframerevcodons), 3):
                            try:
                                parsedic[bc]['codons'].append(inframerevcodons[j:j+3])
                            except KeyError:
                                parsedic[bc]['codons'] = [inframerevcodons[j:j+3]]
                            try:
                                parsedic[bc]['positions'].append((129-((fromend + i[0]+2)//3*3) + j, 129-((fromend + i[0]+2)//3*3) + j + 3))
                            except KeyError:
                                parsedic[bc]['positions'] = [(129-((fromend + i[0]+2)//3*3) + j, 129-((fromend + i[0]+2)//3*3) + j + 3)]
                    fromend += i[0]
    return parsedic



def cigars_to_fasta_errors(cigardic, fastadic, direction):
    bcdic = {}
    for bc in cigardic:
        correctseq = fastadic[cigardic[bc]['name']]
        for cnt, mutants in enumerate(cigardic[bc]['positions']):
            if direction == 'x':
                correctcodon = revcomp(correctseq[mutants[0]:mutants[1]])
                wrongcodon = revcomp(cigardic[bc]['codons'][cnt])
                correctAA = gencode[correctcodon]
                try:
                    wrongAA = gencode[wrongcodon]
                except KeyError:
                    print("barcode ", bc, " with stuff", cigardic[bc], "  is causing issues")

                error = correctAA + str(int((126 - mutants[0])/3)) + wrongAA
            if direction == 'y':
                correctseq = 'G' + correctseq
                correctcodon = correctseq[mutants[0]:mutants[1]]
                wrongcodon = cigardic[bc]['codons'][cnt]
                correctAA = gencode[correctcodon]
                wrongAA = gencode[wrongcodon]
                error = correctAA + str(int(mutants[1] / 3)) + wrongAA
            if bc not in bcdic:
                bcdic[bc] = {}
                bcdic[bc]['name'] = cigardic[bc]['name'] + '__' + error
            else:
                bcdic[bc]['name'] = bcdic[bc]['name'] + '__' + error
        bcdic[bc]['count'] = cnt +1
        bcdic[bc]['unique'] = len(set(cigardic[bc]['positions']))

    return bcdic


def print_bc_dic(bcdic, filename):
    with open(filename, 'w') as f:
        f.write('Barcode,X_peptide,count,unique\n')
        for bc in bcdic:
            mystr = bc + ',' + bcdic[bc]['name'] + ',' + str(bcdic[bc]['count']) + ',' + str(bcdic[bc]['unique']) + '\n'
            f.write(mystr)
    f.close()


if __name__ == '__main__':
    print("Reading x sam")
    cigars_x = samreader('EA66k-full.x.map.sam')
    print("Reading x fasta")
    fastas_x = fastareader('../trimmed_fastas/EA66k-full_x-oligos.fasta', 'x')
    print("parsing x cigar strings")
    pcigs = cigarparse(cigars_x)
    print("turning cigar strings into AAs")
    errors = cigars_to_fasta_errors(pcigs, fastas_x, 'x')
    print("writing x errorfile")
    print_bc_dic(errors, "EA66k.mismatches.x.csv")
    print("doing y's now")
    print("Reading Y sam")
    cigars_y = samreader('EA66k-full.y.map.sam')
    print("Reading y fasta")
    fastay = fastareader('../trimmed_fastas/EA66k-full_y-oligos.fasta', 'y')
    print("parsing cigar strings for y")
    parsycigs = cigarparse_y(cigars_y)
    print("turning cigars strings into AAs for y")
    errorsfory = cigars_to_fasta_errors(parsycigs, fastay, 'y')
    print("Writing y")
    print_bc_dic(errorsfory, 'EA66k.mismatches.y.csv')

    print("pauser")
