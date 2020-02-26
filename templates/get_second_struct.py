import os, re, requests, json
from time import sleep

url = 'http://bioinf.cs.ucl.ac.uk/psipred/api/submission.json'

payload = {'input_data': ('prot.txt', open('prot.txt', 'rb'))}
data = {'job': 'psipred',
        'submission_name': 'testing',
        'email': 'aerijman@fredhutch.org', }

r = requests.post(url, data=data, files=payload)
uid = re.search('{\"UUID\":\"(.*)\",\"sub.*', r.text).group(1)
submission = 'http://bioinf.cs.ucl.ac.uk/psipred/api/submission/'+uid

a = requests.get(submission)
uid = re.search('&quot;/submissions/(.*).mtx&quot;', a.text).group(1)
results = 'http://bioinf.cs.ucl.ac.uk/psipred/api/submissions/'+uid+'.horiz'

not_yet = True
while not_yet:
    
    r = requests.get(results)
    if r.status_code == 200:

        unfiltered = r.text
        break
    sleep(20)   

print(unfiltered)    


#NOTE: Once posted you will need to use the GET submission endpoint
#to retrieve your results. Polling the server about once every 2 or 5 mins
#should be sufficient.
#
# Full details at http://bioinf.cs.ucl.ac.uk/web_servers/web_services/




