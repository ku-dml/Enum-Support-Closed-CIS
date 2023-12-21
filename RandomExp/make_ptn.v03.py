#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 13 16:11:38 2021

@author: student
"""

import random as rand
import sys
import re

def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    return [ atoi(c) for c in re.split(r'(\d+)', text) ]



args = sys.argv
if len(sys.argv)<6:
    sys.stderr.write("python make_ptn.py grh items maxItems outputFile seed\n")
    sys.exit(1)
#grhファイル，q，最大アイテム数，出力ファイル名 seed

graph = args[1]
qlist = range(int(args[2]))
maxI = int(args[3])
output = args[4]
rand.seed(args[5])
file = open(output,mode='w')

ver = []
with open(graph)as f:
    for line in f:
       edge=line.split()
       ver.append(int(edge[0]))
       ver.append(int(edge[2]))

n = list(set(ver))
n.sort()
for i in n:
    file.write(str(i))
    file.write("\t")
    k = rand.randint(0,maxI)
    items = rand.sample(qlist,k)
    if len(items) == 0:
        file.write("!")
    else:
        items.sort()
        for j, e in enumerate(items):
            file.write(str(e))
            if j != len(items)-1:
                file.write(",")
    file.write("\n")
file.close()


