###
# The script will take a series of the multiple genome alignment position information, merge the blocks that are continuous in each species that are within certain range that are defined as "nearby"
# Then the merged blocks can be filtered for certain length cutoff
# By ordering the block accorrding to certain species, the inforamtion of local statics within or across specieis are assigned to an united coordinate system, which is an arbitrary chosen species reference.
###
import sys
#input_pos = sys.argv[1]
#merged_pos = sys.argv[2]
#base = sys.argv[3]
#CLOSE = int(sys.argv[4])
CLOSE = 5000
idx_type = {
    "chr":0,
    "start":1,
    "end":2,
    "direction":3
}
def read_header(header):
    idx_dict = {}
    header_info = header.strip("\n").split("\t")
    for i in range(len(header_info)):
        if header_info[i] != "":
            sp = header_info[i].split("_")[0]
            type = header_info[i].split("_")[1]     
            if sp not in idx_dict:
                idx_dict[sp]=[0,0,0,0]
                idx_dict[sp][idx_type[type]]=i
            else:
                idx_dict[sp][idx_type[type]]=i
    return idx_dict

def continuous_check(chr_current,chr_check,end_current,start_check,direction_current,direction_check):
    if chr_current == chr_check:
        chr_res = 0
    else:
        chr_res = 1
    if abs(start_check - end_current) < CLOSE:
        start_res = 0
    else:
        start_res = 1
    if direction_current == direction_check:
        direction_res = 0
    else:
        direction_res = 1
    if (chr_res +  start_res + direction_res) > 0:
        return 1
    else:
        return 0

def write_current(out,idx_dict,current_chr_dict,current_region_dict,current_direction_dict):
    len_list = []
    for sp in idx_dict:
        start = current_region_dict[sp][0]
        end = current_region_dict[sp][1]
        len_list.append(end - start)
    out.write("{}\t".format(max(len_list)))
    for sp in idx_dict:
        out.write(current_chr_dict[sp]+"\t")
        out.write("{}\t".format(current_region_dict[sp][0]))
        out.write("{}\t".format(current_region_dict[sp][1]))
        out.write(current_direction_dict[sp]+"\t")
    out.write("\n")
    return 0

def init_current(idx_dict,block_info):
    current_chr_dict = {}
    current_region_dict = {}
    current_direct_dict = {}
    for sp in idx_dict:
        current_chr_dict[sp] = block_info[idx_dict[sp][0]]
        current_region_dict[sp] = [int(block_info[idx_dict[sp][1]]),int(block_info[idx_dict[sp][2]])]
        current_direct_dict[sp] = block_info[idx_dict[base][3]]+block_info[idx_dict[sp][3]]
    return current_chr_dict,current_region_dict,current_direct_dict

def update_current(current_chr_dict,current_region_dict,current_direct_dict,idx_dict,block_info):
    for sp in idx_dict:
        current_region_dict[sp][1] = int(block_info[idx_dict[sp][2]])
    return current_chr_dict,current_region_dict,current_direct_dict
    

input_pos = "/home/jilong/spider2/faststorage/social_spiders_2020/people/jilong/scripts/UniCoord/test_data/maf_aligned_pos_demo.tsv"
merged_pos = "merged_demo.tsv"
out = open(merged_pos,"w")
pos_file = open(input_pos)
header_flag = 0
init_flag = 0
base = "DUM"
for line in pos_file:
    if header_flag == 0:
        header = line
        idx_dict = read_header(header)
        header_flag = 1
    else:
        block_info = line.strip("\n").split("\t")
        if init_flag == 0:
            current_chr_dict,current_region_dict,current_direct_dict = init_current(idx_dict,block_info)
            init_flag = 1
        else:
            conti_break = []
            for sp in idx_dict:
                check_chr = block_info[idx_dict[sp][0]]
                check_start = int(block_info[idx_dict[sp][1]])
                check_direction = block_info[idx_dict[base][3]]+block_info[idx_dict[sp][3]]
                current_chr = current_chr_dict[sp]                
                current_end = current_region_dict[sp][1]
                current_direction = current_direct_dict[sp]
                sp_break = continuous_check(current_chr,check_chr,current_end,check_start,current_direction,check_direction)
                conti_break.append(sp_break)
            if sum(conti_break) > 0:            
                write_current(out,idx_dict,current_chr_dict,current_region_dict,current_direct_dict)
                current_chr_dict,current_region_dict,current_direct_dict = init_current(idx_dict,block_info)
            else:
                current_chr_dict,current_region_dict,current_direct_dict = update_current(current_chr_dict,current_region_dict,current_direct_dict,idx_dict,block_info)
