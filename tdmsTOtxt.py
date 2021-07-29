from nptdms import TdmsFile
import glob
import os
import shutil

tdms_file_list=glob.glob("/Users/ahmeduluca/Desktop/Cu/Cu-Raw/**/**/*.tdms")
txt_file_list=glob.glob("/Users/ahmeduluca/Desktop/Cu/Cu-Raw/**/**/*.txt")
process_dir="/Users/ahmeduluca/Desktop/Cu/Cu-Process"
for i in tdms_file_list:
    savefile=os.path.join(process_dir,os.path.basename(os.path.dirname(i)))
    if not os.path.exists(savefile):
        try:
            os.makedirs(savefile)
            print('ok')
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
    for k in txt_file_list:
        if(os.path.basename(os.path.dirname(i)) in os.path.basename(k)):
            shutil.copy(k,savefile)
            cur_dir_txts=glob.glob(os.path.join(os.path.dirname(k)+'/*.txt'))
            break
    with TdmsFile.open(i) as tdms0:
        for group in tdms0.groups():
            savegroup=os.path.join(savefile,group.name)
            print(savegroup)
            if not os.path.exists(savegroup):
                try:
                    os.makedirs(savegroup)
                    print('ok')
                except OSError as exc: # Guard against race condition
                    if exc.errno != errno.EEXIST:
                        raise
            gurup=group.name.replace(" ","")
            try:
                gu, rup=gurup.split('Indentation')
                indent=1
            except:
                if 'Approach' in gurup:
                    gu, rup=gurup.split('Approach')
                    indent=2
                else:
                    gu, rup=gurup.split('Oscillation')
                    indent=0
            for k in cur_dir_txts:
                if(gu in k):
                    if(indent==1 and 'Indentation' in k):
                        shutil.copy(k,savegroup)
                    if(indent==0 and 'Oscillation' in k):
                        shutil.copy(k, savegroup)
                    if(indent==3 and 'Approach' in k):
                        shutil.copy(k,savegroup)
            for channel in group.channels():
                savetxt=os.path.join(savegroup,channel.name+'.txt')
                print(savetxt)
                times=channel.time_track()
                f=open(savetxt,'w+')
                j=0
                for value in channel:
                    f.write(str(times[j]))
                    f.write(" ")
                    f.write(str(value))
                    f.write("\n")
                    j=j+1
                f.close()
                
                          
                          
