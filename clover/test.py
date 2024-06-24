from collections import Counter
from multiprocessing import Process, Queue
import os
import time
from tqdm import tqdm
from clover import load_config as lc
#clover/align module for global alignment not imported

class MyProcess(Process):
    """Process Class for Clover using DUHI's concept"""

    def __init__(self, name, data, q_output):
        Process.__init__(self)
        self.config_dict = lc.out_put_config()
        self.name = name
        self.data = data
        self.q_output = q_output
        self.q_output_temp = []
        self.Cluster_size_threshold = self.config_dict['Cluster_size_threshold'] #unused
        self.ref_list = {}
        self.ref_dict = {}
        self.ref_error_dict = {}
        self.num_dict = {}
        self.tag_dict = {}
        self.index_list = []
        self.now_clust_threshold = self.config_dict['now_clust_threshold']
        self.read_len = self.config_dict['read_len']
        #self.dna_tree_nums = self.config_dict['end_tree_len']
        #self.fuzz_list = [self.config_dict['thd_tree_loc'], self.config_dict['four_tree_loc'], self.config_dict['other_tree_len']]
        self.loc_nums = self.config_dict['Vertical_drift']
        self.tag_nums = self.config_dict['tag_nums']
        self.align_swicth = self.config_dict['align_fuc']
        #self.fuzz_tree_nums = self.config_dict['Horizontal_drift']
        self.h_index = self.config_dict['h_index_nums'] #unused
        self.e_index = self.config_dict['e_index_nums'] #unused
        self.test_num = 0
        self.file_format = "txt"

        # Current implementation only uses .txt input format
        if 'input_path' in self.config_dict: 
            if self.config_dict['input_path'][-1] == "a":
                self.file_format = "fasta"
            elif self.config_dict['input_path'][-1] == "q":
                self.file_format = "fastq"

        # Threshold for inter-cluster distance based on DUHI hash formula - in config file
        self.cluster_hash_threshold = self.config_dict['Cluster_hash_threshold']
        # Initialize arrays for storing cluster hash values
        self.front_hash_array = {}
        self.middle_hash_array = {}
        self.back_hash_array = {}
        self.front_ref_list = {}
        self.middle_ref_list = {}
        self.back_ref_list = {}
    
    def hash_value(self, fragment):
        """DUHI's hash formula for DNA fragment."""
        base_val = {'A': 1, 'T': 2, 'G': 3, 'C': 4}
        hash_val = sum((2 ** (7 - i)) * base_val[char] for i, char in enumerate(fragment))
        return hash_val

    def extract_fragments(self, dna_str):
        """Extract the middle-front, middle, and middle-back fragments of the DNA strand."""
        fragment_length = 8
        middle_start = self.read_len // 3
        middle_end = 2 * middle_start
        
        return (
            dna_str[middle_start:middle_start + fragment_length],
            dna_str[self.read_len // 2 - fragment_length // 2:self.read_len // 2 + fragment_length // 2],
            dna_str[middle_end - fragment_length:middle_end]
        )

    def update_hash(self, dna_num, section):
        hv = 0
        count = 0
        if section == "front":
            for hash_val, hash_count in self.front_ref_list[dna_num]:
                if hash_count > count:
                    hv = hash_val
                    count = hash_count
            self.front_hash_array[dna_num] = (hv, 0)
        elif section == "middle":
            for hash_val, hash_count in self.middle_ref_list[dna_num]:
                if hash_count > count:
                    hv = hash_val
                    count = hash_count
            self.middle_hash_array[dna_num] = (hv, 0)
        elif section == "back":
            for hash_val, hash_count in self.back_ref_list[dna_num]:
                if hash_count > count:
                    hv = hash_val
                    count = hash_count
            self.back_hash_array[dna_num] = (hv, 0)
        else:
            pass

    def cluster(self, read):
        """DUHI's clustering algorithm implemented on Clover."""
        self.test_num += 1
        line_ = read.split()
        
        dna_num = self.test_num
        if self.config_dict['Virtual_mode'] == False or self.config_dict['mmr_mode'] == True:
            dna_index = line_[0]
            dna_str = line_[1]
            dna_tag = line_[0]
        else:
            dna_str = line_[1]
            dna_tag = line_[0]
        #h_index, e_index for dna primer not implemented, assumed 0
        if len(dna_str) < self.config_dict['read_len_min'] or "N" in dna_str:
            return
        
        front, middle, back = self.extract_fragments(dna_str)
        front_hash = self.hash_value(front)
        middle_hash = self.hash_value(middle)
        back_hash = self.hash_value(back)

        # Clustering logic using hash arrays
        # Create new fragments clusters for first dna
        if self.front_hash_array == []:
            self.front_hash_array[dna_num] = (front_hash, 1)
            self.middle_hash_array[dna_num] = (middle_hash, 1)
            self.back_hash_array[dna_num] = (back_hash, 1)
            #self.ref_list[dna_num] = dna_str
            self.front_ref_list[dna_num] = {front_hash: 1}
            self.middle_ref_list[dna_num] = {middle_hash: 1}
            self.back_ref_list[dna_num] = {back_hash: 1}
            self.ref_dict[dna_num] = [dna_tag]
            return
        
        match_found = False
        i = 0
        cluster_num = len(self.front_hash_array)
        for idx in self.front_hash_array:
            count = 0
            front_cluster = self.front_hash_array[idx][0]
            middle_cluster = self.middle_hash_array[idx][0]
            back_cluster = self.back_hash_array[idx][0]
            if front_hash <= front_cluster + self.cluster_hash_threshold and front_hash >= front_cluster - self.cluster_hash_threshold:
                count += 1
            if middle_hash <= middle_cluster + self.cluster_hash_threshold and middle_hash >= middle_cluster - self.cluster_hash_threshold:
                count += 1
            if back_hash <= back_cluster + self.cluster_hash_threshold and back_hash >= back_cluster - self.cluster_hash_threshold:
                count += 1
            
            if count >= 2: # cluster matches if at least 2 fragments meet cluster threshold
                match_found = True
                i = idx
                break

        if match_found:
            # update fragments cluster count
            self.front_hash_array[i][1] += 1
            self.middle_hash_array[i][1] += 1
            self.back_hash_array[i][1] += 1

            # update fragments cluster hash values count
            if front_hash in self.front_ref_list[i]:
                self.front_ref_list[i][front_hash] += 1
            else:
                self.front_ref_list[i][front_hash] = 1
            if middle_hash in self.middle_ref_list[i]:
                self.middle_ref_list[i][middle_hash] += 1
            else:
                self.middle_ref_list[i][middle_hash] = 1
            if back_hash in self.back_ref_list[i]:
                self.back_ref_list[i][back_hash] += 1
            else:
                self.back_ref_list[i][back_hash] = 1

            # update fragment cluster representative hash value if size / 5 = 0
            if self.front_hash_array[i][1] == 5:
                self.update_hash(i, "front")
            if self.middle_hash_array[i][1] == 5:
                self.update_hash(i, "middle")
            if self.back_hash_array[i][1] == 5:
                self.update_hash(i, "back")

            self.ref_dict[i].append(dna_tag)

        else: # if no matches found, create new cluster for all fragments
            self.front_hash_array[dna_num] = (front_hash, 1)
            self.middle_hash_array[dna_num] = (middle_hash, 1)
            self.back_hash_array[dna_num] = (back_hash, 1)
            #self.ref_list[dna_num] = dna_str
            self.front_ref_list[dna_num] = {front_hash: 1}
            self.middle_ref_list[dna_num] = {middle_hash: 1}
            self.back_ref_list[dna_num] = {back_hash: 1}
            self.ref_dict[dna_num] = [dna_tag]
    
    def run(self):
        self.num_dict[self.name + "sum_read_num"] = 0
        self.num_dict[self.name + "error_num"] = 0
        self.num_dict[self.name + "sum_cluster_num"] = 0
        self.num_dict[self.name + "sum_tag"] = []

        if self.config_dict['fast_mode']:
            for line in tqdm(self.data):
                self.cluster(line)
            if 'output_file' in self.config_dict:
                self.num_dict[self.name + "index_list"] = self.index_list

            if self.config_dict['Virtual_mode']:
                tag_sum = 0
                tag_error = 0
                nums_ = 0
                nums_sum = 0
                for key in self.ref_dict:
                    tag_len = len(self.ref_dict[key])
                    if tag_len > self.Cluster_size_threshold:
                        nums_sum += 1
                        tag_c = Counter(self.ref_dict[key]).most_common(1)[0][1]
                        tag_error += tag_len - tag_c
                        tag_sum += tag_len
                        tag = self.ref_dict[key][0]
                        if tag not in self.tag_dict:
                            nums_ += 1
                            self.tag_dict[tag] = 1
                    self.ref_dict[key] = []

                self.num_dict[self.name + "sum_read_num"] += tag_sum
                self.num_dict[self.name + "error_num"] += tag_error
                self.num_dict[self.name + "sum_cluster_num"] = nums_sum
                for key in self.tag_dict:
                    self.num_dict[self.name + "sum_tag"].append(key)
        else:
            if self.config_dict['Virtual_mode']:
                f = open(self.config_dict['input_path'], "r")
                while True:
                    line = f.readline()
                    if line == [] or line == "":
                        break
                    self.cluster(line)

                if 'output_file' in self.config_dict:
                    self.num_dict[self.name + "index_list"] = self.index_list

                tag_sum = 0
                tag_error = 0
                nums_ = 0
                nums_sum = 0
                for key in self.ref_dict:
                    tag_len = len(self.ref_dict[key])
                    if tag_len > self.Cluster_size_threshold:
                        nums_sum += 1
                        tag_c = Counter(self.ref_dict[key]).most_common(1)[0][1]
                        tag_error += tag_len - tag_c
                        tag_sum += tag_len
                        tag = self.ref_dict[key][0]
                        if tag not in self.tag_dict:
                            nums_ += 1
                            self.tag_dict[tag] = 1
                    self.ref_dict[key] = []

                self.num_dict[self.name + "sum_read_num"] += tag_sum
                self.num_dict[self.name + "error_num"] += tag_error
                self.num_dict[self.name + "sum_cluster_num"] = nums_sum
                for key in self.tag_dict:
                    self.num_dict[self.name + "sum_tag"].append(key)
        self.q_output.put(self.num_dict)

def all_permutations(items, length):
    res = []
    def track_back(tmp_permutation):
        if len(tmp_permutation) == length:
            res.append(tmp_permutation[:])
        else:
            for i in items:
                tmp_permutation.append(i)
                track_back(tmp_permutation)
                tmp_permutation.pop()
    track_back([])
    return ["".join(i) for i in res]

def main():
    config_dict = lc.out_put_config()
    tag_nums = config_dict['tag_nums']
    PROCESS_INDEX = int(config_dict['processes_nums'])

    print(config_dict['tag'])

    if PROCESS_INDEX == 0 and config_dict['fast_mode']:
        N_PROCESS = 1
        process_names = ['all']
        print("Generating data")
        data_dict = {'all': []}
        f = open(config_dict['input_path'], 'r').readlines()
        for i in tqdm(f):
            if "N" in i:
                pass
            else:
                data_dict['all'].append(i)
    elif PROCESS_INDEX == 0 and not config_dict['fast_mode']:
        N_PROCESS = 1
        process_names = ['all']
        data_dict = {'all': []}
    elif PROCESS_INDEX > 0:
        N_PROCESS = 4 ** PROCESS_INDEX
        process_names = all_permutations(["A", "T", "G", "C"], PROCESS_INDEX)

        st_pre = time.time()
        data_dict = {}.fromkeys(process_names)
        for i in data_dict:
            data_dict[i] = []

        if config_dict['fast_mode']:
            print("Pre-process the data")
            f = open(config_dict['input_path'], 'r').readlines()
            for i in tqdm(f):
                if "N" in i or "*" in i:
                    pass
                else:
                    data_dict[i.split()[1][:PROCESS_INDEX]].append(i)
        else:
            for i in data_dict:
                data_dict[i] = []

    q_output = Queue(N_PROCESS * 2)

    st = time.time()
    process_dict = {}.fromkeys(process_names)

    for i in process_dict:
        process_dict[i] = MyProcess(i, data_dict[i], q_output)

    for i in process_dict:
        process_dict[i].start()

    st = time.time()
    count_dict = {}
    i = 0
    while True:
        if i == N_PROCESS:
            break
        count_dict.update(q_output.get())
        i += 1
    print("Time:", time.time() - st)

    for i in process_dict:
        process_dict[i].join()

    new_count_dict = {"sum_read_num": 0, "error_num": 0, "sum_cluster_num": 0, "sum_tag": {}, "index_list": []}
    for key in process_dict:
        new_count_dict["sum_read_num"] += count_dict[key + "sum_read_num"]
        new_count_dict["error_num"] += count_dict[key + "error_num"]
        new_count_dict["sum_cluster_num"] += count_dict[key + "sum_cluster_num"]
        tag_list = count_dict[key + "sum_tag"]
        for i in tag_list:
            if i in new_count_dict["sum_tag"]:
                pass
            else:
                new_count_dict["sum_tag"][i] = 1
    if 'output_file' in config_dict and not config_dict["mmr_mode"]:
        for key in process_dict:
            new_count_dict["index_list"] += count_dict[key + "index_list"]
        output_file = open(config_dict['output_file'], 'w')
        output_file.write(str(new_count_dict["index_list"]))
        output_file.close()
    if config_dict['Virtual_mode']:
        print("Number of reads processed:", new_count_dict["sum_read_num"])
        print("Accuracy:", (new_count_dict["sum_read_num"] - new_count_dict["error_num"]) / new_count_dict["sum_read_num"])
        if config_dict['tag_mode']:
            print("Coverage:", len(new_count_dict["sum_tag"]) / int(tag_nums))
            print("Redundancy Rate:", (new_count_dict["sum_cluster_num"] - int(tag_nums)) / int(tag_nums))
    elif config_dict['Statistical_model']:
        pass

    main()
