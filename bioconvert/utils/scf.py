import struct
import copy
from collections import defaultdict

import colorlog
_log = colorlog.getLogger(__name__)


def read_scf(infile):
    sequence = ""
    qualities = []

    with open(infile, "rb") as input_file:

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

    return sequence, qualities, comments



# Return 'length' bits of file 'f_file' starting at offset 'offset'
def read_from_buffer(f_file, length, offset):
    """Return 'length' bits of file 'f_file' starting at offset 'offset'"""
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
                  str(f_file.tell()) + ".")
    return buff


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
    samples = copy.deepcopy(rsamples)
    if direction == "forward":
        for i in range(len(samples)-1, 1, -1):
            samples[i] = samples[i] - 2*samples[i-1] + samples[i-2]
        samples[1] = samples[1] - 2*samples[0]
    elif direction == "backward":
        p_sample1 = 0
        p_sample2 = 0
        for i, _ in enumerate(samples):
            p_sample1 = p_sample1 + samples[i]
            samples[i] = p_sample1 + p_sample2
            p_sample2 = samples[i]
    else:
        msg="Bad direction in 'delta'. Use\" forward\" or\" backward\"."
        _log.critical(msg)
        raise Exception(msg)
    return samples

