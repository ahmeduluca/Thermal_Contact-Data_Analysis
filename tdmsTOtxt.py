from nptdms import TdmsFile
import glob
import os
import shutil

## converts tdms files of experiment to txt & creates new processed file dirs with existing txts
## of exp_log &load &position files & txt converted TDMS group data

# Give the file path to search experiment tdms&txt results:
tdms_file_list=glob.glob("D:/ahmed/RC experiments/Al/**/*.tdms")
txt_file_list=glob.glob("D:/ahmed/RC experiments/Al/**/*.txt")

# Give the new processed file directory path
process_dir="D:/ahmed/RC experiments/Al/Al-Process"

# Walk through found TDMS files:
for i in tdms_file_list:
    
    ## main file dir to save data has date&time of exp on it set via basename:
##! To prevent confusion of days there may be a prior main file for
## each day first -for a  /**/**/.tdms search especially;!- & sorting files for days- 
    savefile=os.path.join(process_dir,os.path.basename(os.path.dirname(i)))
    
    ## Safe creation of new directory if it doesnt exists: 
    if not os.path.exists(savefile):
        try:
            os.makedirs(savefile)
        except OSError as exc:
            if exc.errno != errno.EEXIST:
                raise
    
    ## Copy txt file of exp_logs contained in the experiment files:
    for k in txt_file_list:
        if(os.path.basename(os.path.dirname(i)) in os.path.basename(k)):
            shutil.copy(k,savefile)
            cur_dir_txts=glob.glob(os.path.join(os.path.dirname(k)+'/*.txt'))
            break
        
    ## Use npTDMS' module for reading tdms files: 
    with TdmsFile.open(i) as tdms0:
        
        ## Walk through groups of tdms file that contains steps/calibration of exp:
        for group in tdms0.groups():
            savegroup=os.path.join(savefile,group.name)
            print(savegroup)

            ## Create separate dirs for each step of exp in exp files:
            if not os.path.exists(savegroup):
                try:
                    os.makedirs(savegroup)
                except OSError as exc:
                    if exc.errno != errno.EEXIST:
                        raise

            ## Find appropriate txt files corresponding to Group/step names:
            ## (just to be safe between txt naming and tdms group naming)
            gurup=group.name.replace(" ","")
            try:
                gu, rup=gurup.split('Indentation')
                indent=1
            except:
                if 'Approach' in gurup:
                    gu, rup=gurup.split('Approach')
                    indent=2
                elif 'Retract' in gurup:
                    gu, rup=gurup.split('Retract')
                    indent=4
                elif 'Calibration' in gurup:
                    gu, rup= gurup.split('Calibration')
                    indent=3
                elif 'Oscillation' in gurup:
                    gu, rup=gurup.split('Oscillation')
                    indent=0
                else:
                    gu, rup= "","" ## This TDMS file does not belong to our experiment!
                    print("NOT AN INDENTATION TDMS DATA")
                    break
                
            ## Copy corresponding data found as txt to step directory :
            for k in cur_dir_txts:
                if(gu in k):
                    
##! By directly copying txts, naming for these files left with step name on them with '_'
## separators.. Later it may be standardized for just name of measurement like "Load" etc.
                    if(indent==1 and 'Indentation' in k):
                        shutil.copy(k,savegroup)
                    if(indent==0 and 'Oscillation' in k):
                        shutil.copy(k, savegroup)
                    if(indent==2 and 'Approach' in k):
                        shutil.copy(k,savegroup)
                    if(indent==3 and 'Calibration' in k):
                        shutil.copy(k, savegroup)
                    if(indent==4 and 'Retract' in k):
                        shutil.copy(k, savegroup)

            ## Save tdms data of channels of groups as txt format
            ## with timestamps to step directory:
## (-! while TDMS format logs have low allocation on disk,
## txt converted not. It is just for post process/data visualization 
## with low time consumption for reading..)
            for channel in group.channels():

                ## TXT files have only channel name, bcs each step has its own dir already
                savetxt=os.path.join(savegroup,channel.name+'.txt')
                print(savetxt)

                ## get time track of channel
                times=channel.time_track()
                f=open(savetxt,'w+')
                j=0

                ## write each value of channel log with time track line by line
## ! there should be more time efficient way to do that -with datachunk etc..
                for value in channel:
                    f.write(str(times[j]))
                    f.write(" ")
                    f.write(str(value))
                    f.write("\n")
                    j=j+1
                f.close()
                
                          
                          
