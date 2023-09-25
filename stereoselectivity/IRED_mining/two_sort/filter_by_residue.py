import json


with open("labeled_ireds.json", "r") as f:
    data = json.load(f)

# print(len(data))

# residues = {}
# for i in data:
#     res = i["key_residues"]["187"]
#     if res not in residues.keys():
#         residues[res] = 1
#     else:
#         residues[res] += 1

# print(residues)


###################################### 

# S = []

# c = 0
# for i in data:
#     if i["key_residues"]["139"] == "P":
#         if i["key_residues"]["194"] == "F":
#             if i["key_residues"]["187"] == "Y":
#                 S.append(i)
#                 c += 1

# print(c)

###################################### 

R = []
S = []
U = []

c = 0
for i in data:
    if i["key_residues"]["139"] == "P":
        if i["key_residues"]["194"] == "F":
            if i["key_residues"]["187"] == "Y":
                S.append(i)
            else:
                U.append(i)
        else:
            U.append(i)
    elif i["key_residues"]["139"] == "I" or i["key_residues"]["139"] == "V":
        if i["key_residues"]["194"] == "M" or i["key_residues"]["194"] == "L":
            if i["key_residues"]["187"] == "D":
                R.append(i)
            else:
                U.append(i)
        else:
            U.append(i)
    else:
        U.append(i)

print("S: "+str(len(S)))
print("R: "+str(len(R)))
print("U: "+str(len(U)))

with open("S.fasta", "w+") as f:
    fasta = ""
    for i in S:
        fasta += i["sequence"]
    f.write(fasta)

with open("R.fasta", "w+") as f:
    fasta = ""
    for i in R:
        fasta += i["sequence"]
    f.write(fasta)

with open("U.fasta", "w+") as f:
    fasta = ""
    for i in U:
        fasta += i["sequence"]
    f.write(fasta)