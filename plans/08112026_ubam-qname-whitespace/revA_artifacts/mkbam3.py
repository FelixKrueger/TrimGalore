import struct, sys
CODE = {c:i for i,c in enumerate("=ACMGRSVTWYHKDBN")}
def pack_seq(s):
    out = bytearray()
    for i in range(0, len(s), 2):
        out.append((CODE[s[i]] << 4) | (CODE[s[i+1]] if i+1 < len(s) else 0))
    return bytes(out)
def record(name, seq, flag=4, aux=b""):
    rn = name + b"\x00"
    body = struct.pack("<iiBBHHHiiii", -1,-1,len(rn),0,4680,0,flag,len(seq),-1,-1,0)
    body += rn + pack_seq(seq) + bytes([40])*len(seq) + aux
    return struct.pack("<i", len(body)) + body
seq = "ACGT"*10
hdr = b"@HD\tVN:1.6\n"
out = b"BAM\x01" + struct.pack("<i", len(hdr)) + hdr + struct.pack("<i", 0)
for nm in [b"clean_one", b"bad\nread", b"clean_two"]:
    out += record(nm, seq, 4, b"CBZ" + b"AAACCC\x00")
sys.stdout.buffer.write(out)
