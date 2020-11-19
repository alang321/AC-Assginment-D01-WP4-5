from aircraftProperties import AircraftProperties

temp = []
data = []

file = "data/airfoil" + AircraftProperties.Airfoil["airfoil identifier"] + ".dat"
with open(file) as f:
    lines = f.read().splitlines()

emptyLineIndex = []
#get index of empty lines
for index, line in enumerate(lines):
    if line == "":
        emptyLineIndex.append(index)

temp.append(lines[emptyLineIndex[0]+1:emptyLineIndex[1]])
temp.append(lines[emptyLineIndex[1]+1:])

for index, i in enumerate(temp):
    data.append([[], []])
    for line in i:
        splitline = line.lstrip().split()
        data[index][0].append(float(splitline[0]))
        data[index][1].append(float(splitline[1]))