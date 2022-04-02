import json
import sys
import os


def read_json(json_file):
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data


def outputfile(outputName, content):
    with open(outputName, 'w') as f:
        f.writelines(content)


def transfer_json_to_list(json_file_name, output):
    benchmark_data = read_json(json_file_name)
    challenge_participants = benchmark_data['datalink']['inline_data']['challenge_participants']
    visualization = benchmark_data['datalink']['inline_data']['visualization']
    x_axis, y_axis = visualization['x_axis'], visualization['y_axis']
    evaluate = 'Participant\t%s\t%s\tstderr_x\tstderr_y\n' % (x_axis, y_axis)
    for participant_inf in challenge_participants:
        participant = participant_inf['participant_id']
        x = participant_inf['metric_x']
        y = participant_inf['metric_y']
        x_std = participant_inf['stderr_x']
        y_std = participant_inf['stderr_y']
        evaluate += '%s\t%s\t%s\t%s\t%s\n' % (participant, x, y, x_std, y_std)
    outputfile(output, evaluate)


if __name__ == '__main__':
    fileName = sys.argv[1:][0]
    out = sys.argv[1:][1]
    transfer_json_to_list(fileName, out)