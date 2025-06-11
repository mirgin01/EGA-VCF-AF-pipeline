import csv

def csv_creator(output_csv):
    """
    Create CSV and add proper header
    :param output_csv:
    :return: csv with header and ready to be added new rows
    """
    with open(output_csv, mode="w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["Summary of Samples and Variants"])
        writer.writerow([])  # Empty line for readability


def csv_writer(new_row, output_csv):
    with open(output_csv, mode="a", newline="") as file:
        writer = csv.writer(file)
        for row in new_row:
            writer.writerow(row)