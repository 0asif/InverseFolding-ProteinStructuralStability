import torch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch.nn.functional as F
from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
from get_entropy_probs import get_df_entropies

#NOTE: This is a custom heatmap using stacked plt.bar elements for more control over the coloring
#NOTE: Please use regular heatmaps, this is just to obtain a customized figure

#######################################################################################################
#Functions used in this script
def get_og_counts(df,og_seq):
    count = 0
    retained_og = []
    
    for pos, AA in df.iterrows():
        if og_seq[pos]==AA['top1_AA']:
            count += 1
            retained_og.append(pos)
    return [count, retained_og]

def get_mutants(df,og_seq,l_diff_start):
    l_mutants = []
    
    for pos, AA in df.iterrows():
        mutant = ''
        if og_seq[pos]!=AA['top1_AA']:
            l_mutants.append(og_seq[pos] + str(pos + l_diff_start) + AA['top1_AA'])
    return l_mutants

def get_pos(l):
    l = [int(''.join([char for char in x if char.isdigit()])) for x in l]
    return l

def get_num_matches(known_muts, IF_mutants,label):
    l_IF_mut_pos = [] #store only unique IF mut pos that are close to known mut pos

    for i in get_pos(IF_mutants):
        for j in get_pos(known_muts):
            if ((i-1)==j or (i+1)==j or (i-2)==j or (i+2)==j):
                l_IF_mut_pos.append(i)
    cbp = list(set(l_IF_mut_pos))
    #print('close by pos', len(cbp),'\n')
    
    return cbp

#######################################################################################################
#The actual stuff that gets us the figure
#######################################################################################################

def main():
    #######################################################################################################
    #instantiating data for reference heatmap such as regions on the spike protein, surface exposed residues
    #IUPred dsiorder scores, known mutations etc
    #Also settting important data such as original fasta sequence & relative start position

    infolder = '../data/input/'
    ofolder = '../data/output/'
    l_diff_start = 16 # add this for getting proper mutants with adjusted position

    #reference fasta sequence: EPI_ISL_402124
    fasta = open(infolder + 'spike.fasta','r')
    og_seq = fasta.readlines()[-1].rstrip()
    fasta.close()

    #Run sequence through IUPred and use the 'short' option (https://iupred3.elte.hu/)
    disorder = pd.read_csv(infolder + 'IUPredScores.txt',sep='\t')
    disorder = [x*2 for x in list(disorder['IUPRED SCORE']) if x <= 0.5] #almost all values under 0.5 so we multiply to get better contrast in visuals

    #mutations from: https://www.mdpi.com/article/10.3390/biology13030134/s1 (Table 1)
    mutations = open(infolder + 'spike_mutations.txt','r')
    known_muts = [x.rstrip() for x in mutations.readlines()]
    mutations.close()

    #regions data: https://pmc.ncbi.nlm.nih.gov/articles/PMC7319273/
    wuhan_regions = {
        'NTD':['goldenrod',14,305],
        'N1':['gold',14,26],
        'N2':['gold',67,79],
        'N3':['gold',141,156],
        'N4':['gold',177,186],
        'N5':['gold',246,260],
        'RBD':['darkorange',319,541],
        'RBM':['indianred',437,508],
        'S1/S2':['lavender',685,695],
        'FP':['palevioletred',788,806],
        'HR1':['purple',812,984]
        }

    #Get secondary structure data from Chimera structure viewer
    fh = open(ofolder + 'helices_list.txt','r')
    helices = sorted([int(x) for x in fh.read().rstrip().split(',')])
    fh.close()

    fs = open(ofolder + 'sheets_list.txt','r')
    sheets = sorted([int(x) for x in fs.read().rstrip().split(',')])
    fs.close()

    #list of surface residues: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1010822
    surface = open(infolder + 'exposed_aa.txt','r')
    exposed_aa = [int(x) for x in surface.read().rstrip().split(',')]
    surface.close()

    #Now set the dataframe to get the varying entropy levels
    AA20 = list("SEMDKPTAQRVNGHLIYFWC")
    df = get_df_entropies(ofolder + 'mutation_frequencies_trasposed.csv', AA20)[0]
    
    #This will create different bins for entropy thresholds
    #Minimum entropy goes to 0 and max is under 3
    z01 = df[(df['entropy'] <= 0.01)][['position','top1_AA']]
    z1 = df[(df['entropy'] > 0.01) &(df['entropy'] <= 0.1)][['position','top1_AA']]
    z2 = df[(df['entropy'] > 0.1) & (df['entropy'] <= 0.2)][['position','top1_AA']]
    z3 = df[(df['entropy'] > 0.2) & (df['entropy'] <= 0.3)][['position','top1_AA']]
    z4 = df[(df['entropy'] > 0.3) & (df['entropy'] <= 0.4)][['position','top1_AA']]
    z5 = df[(df['entropy'] > 0.4) & (df['entropy'] <= 0.5)][['position','top1_AA']]
    z10 = df[(df['entropy'] > 0.5) & (df['entropy'] <= 1.0)][['position','top1_AA']]
    z15 = df[(df['entropy'] > 1.0) & (df['entropy'] <= 1.5)][['position','top1_AA']]
    z20 = df[(df['entropy'] > 1.5) & (df['entropy'] <= 2.0)][['position','top1_AA']]
    z25 = df[(df['entropy'] > 2.0) & (df['entropy'] <= 2.5)][['position','top1_AA']]
    z25p = df[df['entropy'] > 2.5]

    dict_pos ={
        '(0.0:0.01]':z01,
        '(0.01:0.1]':z1,
        '(0.1:0.2]':z2,
        '(0.2:0.3]':z3,
        '(0.3:0.4]':z4,
        '(0.4:0.5]':z5,
        '(0.5:1.0]':z10,
        '(1.0:1.5]':z15,
        '(1.5:2.0]':z20,
        '(2.0:2.5]':z25,
        '(2.5:3.0]':z25p
        }

    #Setting up the heatmap figure
    plt.figure(figsize=(14, 8))

    x_pos = [i for i in range(0,1134)]
    y_heights = [1]*1134

    y_i = [y_heights for i in range(len(dict_pos.keys()))]

    y_labels = list(reversed(list(dict_pos.keys())))

    #This section colors the Spike Regions reference bar
    pos_to_color = list(reversed([list(x['position']) for x in dict_pos.values()]))

    bottoms = np.zeros(len(x_pos))

    mut_pos = [(x-1-l_diff_start) for x in list(set(get_pos(known_muts)))]
    for y in range(len(y_i)):

        barlist = plt.bar(x_pos,y_i[y],width=1, bottom=bottoms ,color='white')
        bottoms += y_i[y]

        blue_bars = 0

        for i in pos_to_color[y]:
            barlist[i-1].set_color('blue')
            retained_og = get_og_counts(dict_pos[y_labels[y]], og_seq)
            #-1 is the seconf element of the list as we dont need counts for this one
            #print(retained_og)
            blue_bars += 1
            red_bars = 0
            cyan_bars = 0
        
            for j in retained_og[-1]:
                
                if (j) in mut_pos:
                    barlist[j].set_color('cyan')
                    cyan_bars += 1
                else:
                    barlist[j].set_color('red')
                    red_bars +=1


        plt.text(1200-50,y+0.35,'{:03d}'.format(blue_bars) + ':', color='blue')
        plt.text(1232.5-50,y+0.35,'{:03d}'.format(red_bars) + ':', color='red')
        plt.text(1265-50,y+0.35,'{:03d}'.format(cyan_bars), color='cyan')
                
    plt.text(1145,11,'Total Positions')
    plt.text(1145,10.95,'______________')


    barlist = plt.bar(x_pos,y_i[y],width=1, bottom=bottoms ,color='forestgreen')


    for reg in wuhan_regions.keys():
        temp = wuhan_regions[reg]

        x = (sum(temp[1:])/2)
        y = 11.5
        rotate = 0
        color='white'
        #weight='bold'

        if len(reg)<3 or '/' in reg:
            rotate = 90
            color = 'black'
            x = x-23
            if '/' in reg:
                y=11.15
                x=x+0.6
                
        else:
            x=x-35
        if reg in ['NTD','RBD']:
            x=x-45
        plt.text(x, y, reg, rotation=rotate, color=color)
        
        for pos in range(temp[1],temp[2]):
            barlist[pos-l_diff_start].set_color(temp[0])


    bottoms += y_i[0]

    #This section adds the bar for secondary structural elements
    barlist = plt.bar(x_pos,y_i[0],width=1, bottom=bottoms ,color='bisque')

    for i in helices:
        barlist[i].set_color('rosybrown')

    for i in sheets:
        barlist[i].set_color('maroon')

    #This section will add the bar for surface exposed amino acid positions
    bottoms += y_i[0]

    barlist = plt.bar(x_pos,y_i[0],width=1, bottom=bottoms ,color='darkgoldenrod')

    for i in exposed_aa:
        barlist[i-1].set_color('lightskyblue')


    #Display mutated positions
    bottoms += y_i[0]

    barlist = plt.bar(x_pos,y_i[0],width=1, bottom=bottoms ,color='beige')

    muts = list(set(get_pos(known_muts)))

    for i in muts:
        barlist[i-1-l_diff_start].set_color('black')


    #Display IUPred Scores of positions
    bottoms += y_i[0]

    barlist = plt.bar(x_pos,y_i[0],width=1, bottom=bottoms ,color='black')

    for i in range(len(disorder)):
        barlist[i].set_alpha(disorder[i])

    
    plt.yticks(ticks=[i+0.5 for i in range(16)], labels=y_labels+['Wuhan','α/β/loop','Surface AA','Mutations','Disorder'], fontsize=16)
    xlabels =['16']+[str(x) for x in list(np.arange(50,1101,50))]+['1149']
    xticks = [0]+[i for i in range(34,1101,50)]+[1134]
    plt.xticks(ticks=xticks, labels=xlabels, fontsize=12)

    plt.xlabel("Residue Position", fontsize=20)
    plt.ylabel("Entropy (nats)", fontsize=20)
    plt.title("Original Amino Acids Recovered at Varying Entropy Levels", fontsize=20)



    lw=10
    custom_lines = [Line2D([0], [0], color='lightskyblue', lw=lw, label='Surface AA'),
                    Line2D([0], [0], color='darkgoldenrod', lw=lw, label='Buried'),
                    Line2D([0], [0], color='bisque', lw=lw, label='Loop'),
                    Line2D([0], [0], color='rosybrown', lw=lw, label='α'),
                    Line2D([0], [0], color='maroon', lw=lw, label='β'),
                    Line2D([0], [0], color='red', lw=lw, label='OAA Recovered'),
                    Line2D([0], [0], color='blue', lw=lw, label='OAA Not Recovered'),
                    Line2D([0], [0], color='cyan', lw=lw, label='OAAR but Mutated')]

    plt.legend(handles=custom_lines, loc="upper center", ncol=8).get_frame().set_alpha(0)
    plt.margins(x=0.005)

    plt.savefig('../figures/heatmap.tiff', dpi=150)


if __name__ == "__main__":
    main()