import subprocess

def call_process(command):
    try:
        process=subprocess.call(
            command,
            shell=True
        )
        print(befehl)
        print('ich bin hier Gams')
    except:
        print("Failed {0}".format(process))

if __name__ == '__main__':
    befehl = ['gams','TestLoop.gms','--SORTEDFILE=scenario_sorted_100_1_01']
    call_process(befehl)

