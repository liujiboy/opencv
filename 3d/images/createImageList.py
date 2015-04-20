import os
import sys
root=sys.argv[1]
outfile=open(sys.argv[2],"w")
files=os.listdir(root)
for f in files:
    if f.endswith(".JPG") or f.endswith(".jpg"):
        outfile.write(os.path.join(os.path.abspath(root),f))
        outfile.write("\n")
outfile.close()
