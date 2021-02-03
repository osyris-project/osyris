
import struct


def read_binary_data(content=None, fmt=None, offsets=None, endl=True):
    """
    Unpack binary data from a content buffer using a dict of offsets.
    Also increment the offsets of the corresponding data read, as well as
    increase the line count by 1.
    """

    byte_size = {"b": 1, "h": 2, "i": 4, "q": 8, "f": 4, "d": 8, "e": 8, "n": 8, "l": 8, "s": 1}

    offset = 0
    for key in offsets:
        offset += offsets[key] * byte_size[key]
    # if offset is None:
    #     offset = 4*ninteg + 8*(nlines+nfloat+nlongi) + nstrin + nquadr*16
    if endl:
        offset += 4# + correction

    # byte_size = {"b": 1 , "h": 2, "i": 4, "q": 8, "f": 4, "d": 8, "e": 8}
    if len(fmt) == 1:
        mult = 1
    else:
        mult = int(fmt[:-1])
    pack_size = mult * byte_size[fmt[-1]]

    offsets[fmt[-1]] += mult
    offsets["n"] += 1


    return struct.unpack(fmt, content[offset:offset+pack_size])


def get_binary_data(fmt="",ninteg=0,nlines=0,nfloat=0,nstrin=0,nquadr=0,nlongi=0,offset=None,content=None,correction=0):

    if offset is None:
        offset = 4*ninteg + 8*(nlines+nfloat+nlongi) + nstrin + nquadr*16
    offset += 4 + correction
    byte_size = {"b": 1 , "h": 2, "i": 4, "q": 8, "f": 4, "d": 8, "e": 8}
    if len(fmt) == 1:
        mult = 1
    else:
        mult = eval(fmt[0:len(fmt)-1])
    pack_size = mult*byte_size[fmt[-1]]

    return struct.unpack(fmt, content[offset:offset+pack_size])

