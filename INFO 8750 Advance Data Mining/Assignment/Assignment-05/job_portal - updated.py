import re
from collections import Counter

def main():
    f1 = open('indeed092017.txt','r',errors = 'ignore')
    l = f1.readlines()
    #print(l)
    position =[]
    institution = []
    location = []
    nos = []
    most_job_state = {}
    (nof_job,nof_institute,nof_state) = (0,0,0)
    # pre_line = ''
    #for line in l:
    for i in range(len(l)):
        
        # Search result for rank position
        # (^\w+\s+(Scientist|Professor)) - 20
        # (^[\w+\s+]{0,10}(Scientist|Professor)) - 21
        # (^[\w+\s+]+[\w+\s+]{0,10}(Scientist|Professor)) - 22
        # (^[\w+\s+]{0,10}/{0,1}[\w+\s+]{0,10}/{0,1}(Scientist|Professor)) - 25
        # (^[\w+\s+]{0,10}[or/]{0,1}[\w+\s+]{0,10}[or/]{0,1}[\w+\s+]{0,10}[or/]{0,1}(Scientist|Professor)[\w\s]{0,30}) - 27
        # (^[\w+\s+]{0,10}[,or/]{0,1}[\w+\s+]{0,10}[,or/]{0,1}[\w+\s+]{0,10}[or/]{0,1}(Scientist|Professor)[\s\w\s]{0,30}) -29
        # line 7, 21(save) 66 (1700), 98
        # ---->25

        # Search for Job position
        ra = re.search('(^[\w+\s+]{0,10}[,or/]{0,1}[\w+\s+]{0,10}[,or/]{0,1}[\w+\s+]{0,10}[or/]{0,1}(Scientist|Professor)[\s\w\s]{0,30})', l[i])
        
        # Search result for Institution 
        # (^[\w+\s+]{0,15}[University|Tech][\w+\s+]{0,19}) -29
        # (^[\w+\s+]{0,25}(University|Tech|tech|Hope|UNC)[\w+\s+]{0,19}) - 24
        # (^[\w+\s+&]{0,25}(University|Tech|tech|Hope|UNC)[\w+\s+]{0,19}) - 25

        #ra = re.search('(^[\w\s\.\"]{0,100}(Professor|Director|Coordinator)[\w\s&]{0,40}),',line)

        if ra:
            #print(ra.string)
            #print(ra.group(1))
            #rank.append(ra.group(1))

            # Search result for location
            # [,+\s]{0,3}(\w+\s+)[\d+]{0,3}$ - 25
            # (\w+\s)[,+\s]{0,3}([A-Za-z]{0,10}\s+){0,3}[\d+]{0,6}$ - 25

            # Search for state
            # ([w+\s+\w+]{0,3})[\s+\d+]{0,6}$ - 25
            # Search for location
            # ([\w+\s\w+\s\w+,]{0,20})([w+\s+\w+]{0,3})[\s+\d+]{0,6}$ - 25

            # Search for Institute
            ni = re.search('(^[\w+\s+&]{0,25}(University|Tech|tech|Hope|UNC)[\w+\s+]{0,19})', l[i+1])

            # Search for Location
            #nl = re.search('([\w+\s\w+\s\w+,]{0,20})([w+]{0,20}[\s+\w+]{0,10})[\s+\d+]{0,6}$', l[i+1])
            nl = re.search('([\w+\s\w+\s\w+,]{0,25})([w+]{0,20}[\s+\w+]{0,10})[\s+\d+]{0,6}$', l[i+1])

            # Search for State
            ns = re.search('([\w{0,1}\s{0,1}\w+]{0,3})[\s+\d+]{0,6}$',l[i+1])

            # Task 1: Count the number of Job
            # Check if the real job circular will have tile institute name, and location
            if ra and ni and ns:
                #print(ra.string, ni.string, nl.string, ns.string)
                #print(ns.group(1))
                #print(ra.group(1), ni.group(1), nl.group(1), ns.group(1))

                # Store job posting info
                position.append(ra.group(1))
                institution.append(ni.group(1))
                location.append(nl.group(1))
                nos.append(ns.group(1))
                
                # Counter for job posting
                nof_job+=1


    # Task 1: Print Count the number of Job
    print()
    print("#########################   Task 1: Print Count the number of Job  ###########################")
    print("Total number of Actual Job posting: ", nof_job)
    print()
    
    print("#########################   Task 2: Output the details of the job  ###########################")
    # Task 2: Output the details of the job
    for i in range(nof_job):

        # remove extra newline in the end of each data
        position[i]= position[i].replace("\n", "")
        institution[i]= institution[i].replace("\n", "")
        location[i]= location[i].replace("\n", "")
        nos[i]= nos[i].replace("\n", "")
        print ("Job Post #"+ str(i+1)+" : "+position[i]+", "+institution[i]+ ", "+location[i])
    
    print()
    # Task 3: States have the most jobs
    print("#########################   Task 3:  States have the most jobs  ###########################")
    most_job_state = dict(Counter(nos))

    print("State wise Job posting (Unsorted):\n", most_job_state)
    print()

    #max_key = max(most_job_state, keys = most_job_state.get)
    max_key = [k for k, v in most_job_state.items() if v == max(most_job_state.values())]
    max_value = max(most_job_state.values())
    print("Find which state(s) have the most jobs: ")
    print(max_key," has "+str(max_value)+" jobs in total.")
    print()

    # Descending Sorting based on number of state job posting 
    sorted_keys = sorted(most_job_state, key=most_job_state.get, reverse=True)
    print("State wise most jobs(Sorted):")
    for r in sorted_keys:
        print(r," - ", most_job_state[r])









main()