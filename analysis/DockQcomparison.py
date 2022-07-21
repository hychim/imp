#!/usr/bin/env python3

import pandas as pd
import sys
import os
import argparse
import subprocess
import re
import numpy as np
from Bio.PDB import PDBIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBParser import PDBParser
import string


def getname(id1,id2):
    if (id1>id2):
        Name=id2+"-"+id1
    else:
        Name=id1+"-"+id2
    return (Name)

def structtonumpy(structure):
    chains={}
    chain_ids=[]
    for model in structure:
        for chain in model:
            ch=chain.get_id()
            chain_ids+=[ch]
            chains[ch] = []
            for residue in chain:
                try:
                    x,y,z=residue['CB'].get_vector()
                    chains[ch].append([x,y,z])
                except:
                    try:
                        x,y,z=residue['CA'].get_vector()
                        chains[ch].append([x,y,z])
                    except:
                        try:
                            x,y,z=residue['N'].get_vector()
                            chains[ch].append([x,y,z])
                        except:
                            i=0
    return (chains,chain_ids)



parser = argparse.ArgumentParser(description = '''Calculate DockQ valued for corrsponding pairsin two complexes.''')

parser.add_argument("-p","--pdb",type=argparse.FileType('r'),help="Input pdb file",required=True)
parser.add_argument("-m","--model",type=argparse.FileType('r'),help="Input model file",required=True)
parser.add_argument("-a","--assembly",type=argparse.FileType('r'),help="Input Assembly path",required=False)
parser.add_argument("-C","--cutoff", type=float,required=False,default=8.0,help="Cutoff for defining distances")
parser.add_argument("-min","--minnum", type=float,required=False,default=20.0,help="Minimum Number of contacts to be included")
parser.add_argument("-n","--nommalign", action='store_true',help="skip MMalign to match chains")
parser.add_argument("-M","--mmalign", type=str,required=False,default="/home/arnee/git/huintaf2/bin/MMalign",help="Path to MMalign")
parser.add_argument("-D","--DockQ", type=str,required=False,default="/home/arnee/git/DockQ/DockQ.py",help="Path to DockQ")
parser.add_argument("-o","--output", type=str,required=False,default="",help="Output CSV File")





args = parser.parse_args()

MM=args.mmalign
DockQ=args.DockQ



#  First we need to match chains (for homologous sequences)
matching={}
if (args.nommalign==False):
    MMout = subprocess.check_output([MM,args.pdb.name, args.model.name])
    MMout = MMout.decode('utf8').split("\n")
    for line in MMout:
        #print (line)
        if line.find("Name of Structure_1:")!=-1:
            line= line.replace(' (to be superimposed onto Structure_2)','')
            chains1=line.split(":")
            #chains1.pop
            #chains1.pop

        elif line.find("Name of Structure_2:")!=-1:
            chains2=line.split(":")
            #chains2.pop
            #chains2.pop


    for i in range(len(chains1)):
        matching[chains1[i]]=chains2[i]
else:
    for i in string.printable:
        matching[i]=i

#print (chains1)
#print (chains2)
#print (matching)

# secondly we need to find all "contacting chains in original PDB"
Name=args.pdb.name
bio_parser = PDBParser()
structure_file = args.pdb
structure_id = args.pdb.name[:-4][-4:]
structure = bio_parser.get_structure(structure_id, structure_file)

cutoff=args.cutoff
minnum=args.minnum
CHAINS,IDS=structtonumpy(structure)
#m = np.sqrt(np.sum((CA[:,np.newaxis,:] -CA[np.newaxis,:,:])**2, axis=2))
#n=np.where(m<=cutoff)
#print (len(CHAINS),len(IDS))
#print (CHAINS,IDS)
#print ("Ch1,Ch2,NumContacts")
df_results=pd.DataFrame(columns=['PDBchain1', 'PDBchain2', 'Modelchain1',"Modelchain2","Numcontacts","DockQ"])
for i in IDS:
    for j in IDS:
        if ( j > i):
            #Get interface
            coords1, coords2 = np.array(CHAINS[i]), np.array(CHAINS[j])
            #Check total length
            #Calc 2-norm
            mat = np.append(coords1, coords2,axis=0)
            a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
            dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T
            l1 = len(coords1)
            contact_dists = dists[:l1,l1:]
            contacts = np.argwhere(contact_dists<=cutoff)
            #print (str(i)+","+str(j)+","+str(len(contacts)))
            if (len(contacts)>minnum):
                try:
                    cmd="/usr/bin/python3 "+DockQ +" "+ args.pdb.name +" "+ args.model.name + " -model_chain1 " + matching[i] + " -model_chain2 " + matching[j] + " -native_chain1 " + i + " -native_chain2 " +j
                    DockQout = subprocess.check_output(cmd.split())
                    DockQout = DockQout.decode('utf8').split("\n")
                    #print (DockQout.stdout)
                    #print (DockQout.stderr)
                    for line in DockQout:
                        #print (line)
                        if line.startswith("DockQ"):
                            #print (line)
                            x,y=line.split()
                    print (i,j,len(contacts),matching[i],matching[j],y)
                    df_results=df_results.append({'PDBchain1':str(i), 'PDBchain2':str(j), 'Modelchain1':str(matching[i]),"Modelchain2":str(matching[j]),"Numcontacts":len(contacts),"DockQ":y}, ignore_index=True)
                except:
                    print('Extra interface')

# Read optimal path if available
df_results["PDB"]=structure_id

if args.assembly:
    df_assembly=pd.read_csv(args.assembly.name)
    df_assembly["Name"] = df_assembly.apply(lambda x: getname(x["Chain"],x["Edge_chain"]),axis=1)
    df_results["Name"] = df_results.apply(lambda x: getname(x["Modelchain1"],x["Modelchain2"]),axis=1)
    df_results=pd.merge(df_assembly,df_results,on="Name",how="outer")

if (args.output!=""):
    df_results.to_csv(args.output,index=False)
