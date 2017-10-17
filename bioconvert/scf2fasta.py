import sys
import struct
import copy
from collections import defaultdict
from bioconvert import ConvBase

class Scf2Fasta(ConvBase):
    """
    Converts a binary scf file to Fasta format.

    :param str infile:
    :param str outfile:
    """
    input_ext = '.scf'
    output_ext = ['.fa', '.fasta']

    def __call__(self, *args, **kwargs):
        sequence = ""
        qualities = []

        with open(self.infile, "rb") as input_file:

            # Header
            """
            char[4] magic_number -> usually .scf
            uint samples -> Number of elements in Samples matrix
            uint samples_offset -> Byte offset from start of file
            uint bases -> Number of bases in Bases matrix
            uint bases_left_clip -> OBSOLETE: No. bases in left clip (vector)
            uint bases_right_clip -> OBSOLETE: No. bases in right clip (qual)
            uint bases_offset -> Byte offset from start of file
            uint comments_size -> Number of bytes in Comment section
            uint comments_offset -> Byte offset from start of file
            char[4] version -> "ver.rev", eg '2' '.' '0' '0' or '3' '.' '0' '0'
            uint sample_size -> Size of samples in bytes 1=8bits, 2=16bits*/
            uint code_set -> code set used (but ignored!)*/
            uint private_size -> No. of bytes of Private data, 0 if none
            uint private_offset -> Byte offset from start of file
            char[18] spare -> Unused
            """
            buff = input_file.read(56)
            magic_number, samples, samples_offset, bases, bases_left_clip, \
            bases_right_clip, bases_offset, comments_size, comments_offset, \
            version, sample_size, code_set, private_size, \
            private_offset = struct.unpack("!4s 8I 4s 4I", buff)
            # Header is supposed to be 128B
            if samples_offset != 128:
                print("Warning: possible bad file")
            # End of header, unused data
            spare_size = samples_offset - 56 # Should be 72
            buff = input_file.read(spare_size)
            #spare = struct.unpack("!%dI" % (spare_size/4), buff)

            # Trying to replicate the next_seq perl function
            # V2
            if float(version) < 3:

                accuracies = defaultdict(list)
                qualities = []
                length = samples * sample_size * 4
                buff = read_from_buffer(input_file, length, samples_offset)

                # Traces need to be split
                #traces = parse_v2_traces(buff, sample_size)

                # Get the base information
                length = bases * 12
                input_file.seek(bases_offset)
                buff = read_from_buffer(input_file, length, bases_offset)
                # Distill the information into its fractions
                for i in range(0, length, 12):
                    # uint peak_index -> Index into Samples matrix for base posn
                    # uchar prob_A -> Probability of it being an A
                    # uchar prob_C -> Probability of it being an C
                    # uchar prob_G -> Probability of it being an G
                    # uchar prob_T -> Probability of it being an T
                    # char base -> Called base character
                    # uchar[3] spare -> Spare
                    read = struct.unpack("!I B B B B s 3B", buff[i:i+12])
                    #indices = read[0]
                    curbase = read[5].lower().decode("utf-8")
                    if curbase == "a":
                        currqual = read[1]
                        accuracies[curbase].append(currqual)
                    elif curbase == "c":
                        currqual = read[2]
                        accuracies[curbase].append(currqual)
                    elif curbase == "g":
                        currqual = read[3]
                        accuracies[curbase].append(currqual)
                    elif curbase == "t":
                        currqual = read[4]
                        accuracies[curbase].append(currqual)
                    else:
                        currqual = -1

                    sequence += str(curbase).upper()
                    qualities.append(currqual)

            # V3
            else:
                length = samples * sample_size
                transformed_read = []
                # Just loop 4 times
                for _ in "acgt":
                    buff = read_from_buffer(input_file, length, samples_offset)
                    byte = "!%dH" % samples
                    # Not tested
                    if sample_size == 1:
                        byte = "%db" % samples

                    read_tmp = struct.unpack(byte, buff)

                    # this little spurt of nonsense is because
                    # the trace values are given in the binary
                    # file as unsigned shorts but they really
                    # are signed deltas. 30000 is an arbitrary number
                    # (will there be any traces with a given
                    # point greater then 30000? I hope not.
                    # once the read is read, it must be changed
                    # from relative
                    read = []
                    for i in read_tmp:
                        if i > 30000:
                            i -= 65536
                        read.append(i)
                    transformed_read += delta(read, "backward")

                    # For 8-bit data we need to emulate a signed/unsigned
                    # cast that is implicit in the C implementations
                    # Not tested
                    if sample_size == 1:
                        for i, _ in enumerate(transformed_read):
                            if transformed_read[i] < 0:
                                transformed_read[i] += 256

                # Get the peak index information
                offset = bases_offset
                length = bases * 4
                buff = read_from_buffer(input_file, length, offset)
                #peak_indices = get_v3_peak_indices(buff)
                offset += length
                # Get the accuracy information
                buff = read_from_buffer(input_file, length, offset)
                accuracies = get_v3_base_accuracies(buff)
                # Get the base information
                offset += length
                length = bases
                buff = read_from_buffer(input_file, length, offset)
                seq = struct.unpack("%ds" % length, buff)
                for i in seq:
                    sequence += i.decode("utf-8")
                # Extract the calls from the accuracy information.
                qualities = get_v3_quality(sequence, accuracies)
                # Bug in V3, 1 bytes is added for unknown reason
                comments_size = comments_size - 1

            # Comments
            comments = ""
            input_file.seek(comments_offset)
            # 1 bytes is added to comments for unknown reason (2 in V3)
            while input_file.tell() < comments_offset + comments_size - 1:
                tmp = input_file.read(1)
                comments += struct.unpack('c', tmp)[0].decode("utf-8")

        # Wrinting output file
        with open(self.outfile, "w") as output_file:
            output_file.write(">" + comments.replace("\n", "-").replace(" ", "_") + "\n")
            output_file.write(sequence)

        """
        print(sequence)
        print("INFORMATIONS")
        print("magic_number = " + magic_number.decode("utf-8"))
        print("samples = " + str(samples))
        print("samples_offset = " + str(samples_offset))
        print("bases = " + str(bases))
        print("bases_left_clip = " + str(bases_left_clip))
        print("bases_right_clip = " + str(bases_right_clip))
        print("bases_offset = " + str(bases_offset))
        print("comments_size = " + str(comments_size))
        print("comments_offset = " + str(comments_offset))
        print("version = " + version.decode("utf-8"))
        print("sample_size = " + str(sample_size))
        print("code_set = " + str(code_set))
        print("private_size = " + str(private_size))
        print("private_offset = " + str(private_offset))
        #print("spare = " + str(spare))
        print()
        print("COMMENTS")
        print(comments)
        print()
        print("SEQUENCE")
        print(sequence)
        print()
        print("QUALITIES")
        print(qualities)
        """




"""
http://staden.sourceforge.net/manual/formats_unix_2.html
http://doc.bioperl.org/bioperl-live/Bio/SeqIO/scf.html#POD6
https://docs.python.org/2/library/struct.html
http://www.perlmonks.org/?node_id=224666

SCF file organisation (more or less)

Length in bytes                        Data
---------------------------------------------------------------------------
128                                    header
Number of samples * sample size        Samples for A trace
Number of samples * sample size        Samples for C trace
Number of samples * sample size        Samples for G trace
Number of samples * sample size        Samples for T trace
Number of bases * 4                    Offset into peak index for each base
Number of bases                        Accuracy estimate bases being 'A'
Number of bases                        Accuracy estimate bases being 'C'
Number of bases                        Accuracy estimate bases being 'G'
Number of bases                        Accuracy estimate bases being 'T'
Number of bases                        The called bases
Number of bases * 3                    Reserved for future use
Comments size                          Comments
Private data size                      Private data
"""

# Return 'length' bits of file 'f_file' starting at offset 'offset'
def read_from_buffer(f_file, length, offset):
    f_file.seek(offset)
    buff = f_file.read(length)
    if len(buff) != length:
        print("The read was incomplete! Trying harder.")
        missing_length = length - len(buff)
        buff2 = f_file.read(missing_length)
        buff += buff2
        if len(buff) != length:
            print("Unexpected end of file while reading from SCF file. I \
                  should have read " + str(length) + " but instead got " +
                  str(len(buff)) + "! Current file position is " +
                  f_file.tell() + ".")
    return buff

# Parses a v2 scf trace array into its base components.
def parse_v2_traces(buff, sample_size):
    traces = []
    # Not tested
    byte = "!%dH" % len(buff)
    # Tested
    if sample_size == 1:
        byte = "%db" % len(buff)

    read = struct.unpack(byte, buff)
    # Traces for a, c, g and t
    traces = defaultdict(list)
    for i in range(0, len(read), 4):
        traces['a'].append(read[i])
        traces['c'].append(read[i+1])
        traces['g'].append(read[i+3])
        traces['t'].append(read[i+2])
    return traces

# Unpacks the peak indices accuracies for v3 scf
def get_v3_peak_indices(buff):
    length = len(buff)
    read = struct.unpack('!%dI' % (length/4), buff)[0]
    return read

# Unpack the base accuracies for v3 scf
def get_v3_base_accuracies(buff):
    # 4 keys: a, c, g and t
    accuracies = defaultdict(list)
    length = len(buff)
    qlength = int(length/4)
    offset = 0
    # For all bases
    for curbase in "acgt":
        read = struct.unpack('%dB' % qlength, buff[offset:offset+qlength])
        for i in read:
            accuracies[curbase].append(i)
        offset += qlength
    return accuracies

# Get the quality for all bases in v3 scf
def get_v3_quality(sequence, accuracies):
    qualities = []
    for i, _ in enumerate(sequence):
        curbase = sequence[i].lower()
        if curbase == "a" or curbase == "c" or curbase == "g" or curbase == "t":
            currqual = accuracies[curbase][i]
        else:
            # N value most of the time
            currqual = -1
        qualities.append(currqual)
    return qualities

# If job == DELTA_IT:
#     Change a series of sample points to a series of delta delta values:
#     ie change them in two steps:
#     first: delta = current_value - previous_value
#     then: delta_delta = delta - previous_delta
# else
#     do the reverse
def delta(rsamples, direction):

    slow_but_clear = 0

    samples = copy.deepcopy(rsamples)

    if direction == "forward":
        if slow_but_clear:
            p_delta = 0
            for i, _ in enumerate(samples):
                p_sample = samples[i]
                samples[i] = samples[i] - p_delta
                p_delta = p_sample
            p_delta = 0
            for i, _ in enumerate(samples):
                p_sample = samples[i]
                samples[i] = samples[i] - p_delta
                p_delta = p_sample
        else:
            for i in range(len(samples)-1, 1, -1):
                samples[i] = samples[i] - 2*samples[i-1] + samples[i-2]

            samples[1] = samples[1] - 2*samples[0]

    elif direction == "backward":
        if slow_but_clear:
            p_sample = 0
            for i, _ in enumerate(samples):
                samples[i] = samples[i] + p_sample
                p_sample = samples[i]
            p_sample = 0
            for i, _ in enumerate(samples):
                samples[i] = samples[i] + p_sample
                p_sample = samples[i]
        else:
            p_sample1 = 0
            p_sample2 = 0
            for i, _ in enumerate(samples):
                p_sample1 = p_sample1 + samples[i]
                samples[i] = p_sample1 + p_sample2
                p_sample2 = samples[i]

    else:
        print("Bad direction in 'delta'. Use\" forward\" or\" backward\".")
        sys.exit(1)
    return samples
