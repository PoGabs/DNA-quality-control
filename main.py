from collections import Counter
import re
import gzip
# get sum of negative aspects, the lowest sum wins
file1 = gzip.open(input())
file2 = gzip.open(input())
file3 = gzip.open(input())
list_of_lists = []
app_counter = 0
while app_counter < 3:
    user_input = ""
    before_read = "@SRR"
    if app_counter == 0:
        user_input = file1
    elif app_counter == 1:
        user_input = file2
    elif app_counter == 2:
        user_input = file3
    file_content = (user_input.read()).decode("utf-8")
    file_content_lines = file_content.split("\n")
    found = []
    n = 0
    for i in file_content_lines:
        if i[:4] == before_read:
            found.append(n)
        n = n + 1


    def get_read(found_reads):
        read_lines = []
        for read in found_reads:
            read_lines.append(read + 1)
        return read_lines


    def get_read_lines(reads):
        line_set = set()
        for a in reads:
            line_set.add(file_content_lines[a].rstrip("\n"))
        return list(line_set)


    def determine_lengths(reads):
        total_number_of_reads = len(reads)
        line_list = get_read_lines(reads)
        lens = []
        for b in line_list:
            lens.append(len(b))
        lens_frequency = Counter(lens)
        sum_of_lengths = 0
        for d in sorted(lens_frequency):
            sum_of_lengths = sum_of_lengths + (d * lens_frequency[d])
        sum_of_occurrences = 0
        for e in sorted(lens_frequency):
            sum_of_occurrences = sum_of_occurrences + lens_frequency[e]
        average = round(sum_of_lengths / sum_of_occurrences)
        return [total_number_of_reads, average]


    def g_c_percentage(reads, content_lines):
        read_lines = []
        for number in reads:
            read_lines.append(content_lines[number])
        average_list = []
        for a in read_lines:
            g_count = 0
            c_count = 0
            a_count = 0
            t_count = 0
            n_count = 0
            for nucleotide in a:
                if nucleotide == "G":
                    g_count = g_count + 1
                elif nucleotide == "C":
                    c_count = c_count + 1
                elif nucleotide == "A":
                    a_count = a_count + 1
                elif nucleotide == "T":
                    t_count = t_count + 1
                elif nucleotide == "N":
                    n_count = n_count + 1
            above_line = g_count + c_count
            below_line = a_count + t_count + g_count + c_count + n_count
            above_divided_below = above_line / below_line
            average_list.append(above_divided_below)
        complete_mean = sum(average_list) / len(average_list)
        return round(complete_mean * 100, 2)


    def find_repeat1(reads, file_lines):
        full_reads = []
        for number in reads:
            full_reads.append(file_lines[number])
        result = Counter(full_reads)
        repeating = []
        for value in result:
            if result[value] > 1:
                repeating.append(result[value])
        calculation = sum(repeating) - len(repeating)
        return calculation


    def counts_with_n(reads, file_lines):
        full_reads = []
        for number in reads:
            full_reads.append(file_lines[number])
        n_counter = 0
        for read in full_reads:
            if bool(re.search("N", read)):
                n_counter = n_counter + 1
        return n_counter


    def n_mean(reads, file_lines):
        full_reads = []
        for number in reads:
            full_reads.append(file_lines[number])
        average_list = []
        for a in full_reads:
            n_count = 0
            for nucleotide in a:
                if nucleotide == "N":
                    n_count = n_count + 1
            above_line = n_count
            below_line = len(a)
            result = 100 * (above_line / below_line)
            average_list.append(result)
        complete_mean = sum(average_list) / len(average_list)
        return round(complete_mean, 2)


    given_reads = get_read(found)
    main_read_data = determine_lengths(given_reads)
    reads_in_file = main_read_data[0]
    reads_sequence_average_length = main_read_data[1]
    repeat_number = find_repeat1(given_reads, file_content_lines)
    reads_with_n = counts_with_n(given_reads, file_content_lines)
    g_c_content = g_c_percentage(given_reads, file_content_lines)
    n_per_read_sequence = n_mean(given_reads, file_content_lines)
    user_input.close()
    list_of_lists.append([reads_in_file, reads_sequence_average_length, repeat_number, reads_with_n, g_c_content, n_per_read_sequence])
    app_counter = app_counter + 1


def quality_comparison(lists_list):
    first_negatives = [lists_list[0][3], lists_list[0][5]]
    second_negatives = [lists_list[1][3], lists_list[1][5]]
    third_negatives = [lists_list[2][3], lists_list[2][5]]
    list_for_comparison = [sum(first_negatives), sum(second_negatives), sum(third_negatives)]
    best = min(list_for_comparison)
    good_index = list_for_comparison.index(best)
    return good_index


best_read_no = quality_comparison(list_of_lists)
print("Reads in the file = " + str(list_of_lists[best_read_no][0]))
print("Reads sequence average length = " + str(list_of_lists[best_read_no][1]))
print()
print("Repeats = " + str(list_of_lists[best_read_no][2]))
print("Reads with Ns = " + str(list_of_lists[best_read_no][3]))
print()
print("GC content average = " + str(list_of_lists[best_read_no][4]) + "%")
print("Ns per read sequence = " + str(list_of_lists[best_read_no][5]) + "%")
