# This script converts the commands in polycomb_var.txt into a series of bash
# scripts which upload the results to cloud storage.
with open("PolycombCode/src/polycomb_var.txt", "r") as f:
    commands = f.readlines()

for command in commands:
    script = []
    tokens = command.split(" ")
    number = tokens[-1] # This is 2020_1, 2020_2, etc.
    number = number.strip()
    script.append("cd /polycomb\n")
    script.append("num={0}\n".format(number))
    script.append("echo \"Running job for $num\"\n")
    script.append("{0}".format(command))
    script.append("python3 upload.py polycomb $num results/$num\n")
    with open("PolycombCode/scripts/{0}.sh".format(number), "w") as f:
        f.writelines(script)
