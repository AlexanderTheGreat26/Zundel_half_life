import re

bad_words = ['#\n', " 2 1 ", " 3 2 ", " 1 0 ", " 2 0 ", " 3 0 ", " 1 1 ", " 1 2 ",
             "# Number of particles", "# Max number of bonds per atom ",
             "# Particle connection table and bond orders",
             "# id type nb id_1...id_nb mol bo_1...bo_nb abo nlp q"]

with open ('bonds.reaxff', 'r') as f, open("buf", 'w') as buf:
    for line in f:
        if not any(bad_word in line for bad_word in bad_words):
            buf.write(line)

with open("buf", 'r') as buf:

    content = buf.read()
    
    content_new = re.sub("         *.*", "", content, flags = re.M)

    content_new = re.sub(" *\n", "\n", content_new, flags = re.M);

with open("bonds.only", "wt", encoding="utf-8") as file_handle:
    file_handle.write(content_new)
