import csv
csv1 = "C:\\Users\\Daradis\\Desktop\\test\\test1.csv"
csv2 = "C:\\Users\\Daradis\\Desktop\\test\\test2.csv"
commons_from_csv1 = "C:\\Users\\Daradis\\Desktop\\test\\commons_from_csv1.csv"
commons_from_csv2 = "C:\\Users\\Daradis\\Desktop\\test\\commons_from_csv2.csv"
DE_same_direction_from_csv1 = "C:\\Users\\Daradis\\Desktop\\test\\DE_same_direction_from_csv1.csv"
DE_same_direction_from_csv2 = "C:\\Users\\Daradis\\Desktop\\test\\DE_same_direction_from_csv2.csv"
ID_COLUMN_INDEX = 0  # some ID (ensembl, HUGO)
VALUE_COLUMN_INDEX = 3  # fold expression column


#This defines an argument "commons" which is a list containing shared items between csv1 and csv2 in row specified by ID_COLUMN_INDEX
def read_ids_from_file(file_path):
    with open(file_path, "r") as file:
        reader = csv.reader(file)
        next(reader)
        return [row[ID_COLUMN_INDEX] for row in reader]


#This function creates a csv file which doesn't contain item that are not shared between csv1 and csv2
def write_common_rows(source_file, commons_from_csv, common_ids):
    with open(source_file, "r") as source, open(commons_from_csv, "w", newline="") as result:
        reader = csv.reader(source)
        writer = csv.writer(result)
        next(reader)
        writer.writerows(row for row in reader if row[0] in common_ids)

#This function creates a list of dictionaries. Key specified by ID_COLUMN_INDEX
#The value must be a numeric value that will have its sign compared – VALUE_COLUMN_INDEX.
def read_dict_from_file(file_path):
    with open(file_path, "r") as file:
        reader = csv.reader(file)
        return [{row[ID_COLUMN_INDEX]: row[VALUE_COLUMN_INDEX]} for row in reader]

#This function compares and filters out values in column specified by function reads_ids_from_file
#requirements to pass the filter are: that numberical values in both csv has the same sign; second req
#checks if the ID of the to be written row corresponds to the value
def write_rows_matching_condition(source_file, DE_same_direction_from_csv):
    with open(source_file, "r") as source, open(DE_same_direction_from_csv, "w", newline="") as result:
        reader = csv.reader(source)
        writer = csv.writer(result)
        
        for i, j, row in zip(file1_ids, file2_ids, reader):
            if float(list(i.values())[0]) > 0 and float(list(j.values())[0]) > 0 and row[ID_COLUMN_INDEX] == list(i.keys())[0]:
                writer.writerow(row)
            elif float(list(i.values())[0]) < 0 and float(list(j.values())[0]) < 0 and row[ID_COLUMN_INDEX] == list(i.keys())[0]:
                writer.writerow(row)


#####
file1_ids = read_ids_from_file(csv1)
file2_ids = read_ids_from_file(csv2)
commons = set(file1_ids).intersection(file2_ids)
#####
#####
write_common_rows(csv1, commons_from_csv1, commons)
write_common_rows(csv2, commons_from_csv2, commons)
#####
#####
file1_ids = read_dict_from_file(commons_from_csv1)
file2_ids = read_dict_from_file(commons_from_csv2)
#####
#####
write_rows_matching_condition(commons_from_csv1, DE_same_direction_from_csv1)
write_rows_matching_condition(commons_from_csv2, DE_same_direction_from_csv2)
#####



"""             read_ids_from_file
    Reads a specific column (ID column) from a CSV file and returns the IDs as a list.

    Parameters:
        file_path (str): Path to the CSV file to read from.

    Returns:
        A list containing the IDs from the specified column.

    Behavior:
        The function reads the contents of the CSV file specified by file_path.
        It assumes that the first row of the file contains headers and skips it.
        The function extracts the values from the specified ID column (ID_COLUMN_INDEX)
        and returns them as a list.
"""


"""                 write_common_rows
    Writes rows from the source_file to the commons_from_csv file,
    filtering out rows whose first element (ID) is not present in the common_ids set.

    Parameters:
        source_file (str): Path to the source CSV file to read from.
        commons_from_csv1 (str): Path to the target CSV file to write the filtered rows.
        common_ids (set): A set containing the common IDs shared between two CSV files.

    Behavior:
        The function reads rows from the source_file and writes the rows to the commons_from_csv file
        if their first element (ID) is present in the common_ids set.
"""


"""               read_dict_from_file
    Reads two specific columns from a CSV file and returns the data as a list of dictionaries.

    Parameters:
        file_path (str): Path to the CSV file to read from; usually file from write_common_rows.

    Returns:
        list: A list containing dictionaries with data from specified columns.

    Behavior:
        The function reads the contents of the CSV file specified by file_path.
        The function creates a dictionary for each row, where the key is extracted from the
        specified ID column (ID_COLUMN_INDEX) and the value from the specified VALUE column (VALUE_COLUMN_INDEX).
        The dictionaries are collected in a list and returned.

"""


"""                  write_rows_matching_condition         

    Writes rows from the source_file to the commons_from_csv1 file,
    based on specific conditions and comparisons with dictionaries
    from the two common ID files.

    Parameters:
        source_file (str): Path to the source CSV file to read from.
        DE_same_direction_from_csv (str): Path to the target CSV file to write the filtered rows.

    Behavior:
        The function reads rows from the source_file and compares values with dictionaries
        in file1_ids and file2_ids. It writes the rows to the commons_from_csv1 file if the
        following conditions are met:
        1. The numeric values in both dictionaries have the same sign (positive or negative).
        2. The ID in the row matches the key in the dictionary.
"""