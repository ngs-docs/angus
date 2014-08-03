1. Edit conf.py and replace 'labibi' in the 'project' line.

2. Either replace or eliminate the Google Analytics ID, the disqus name,
   and the github information. ::

     html_context = {
       "google_analytics_id" : 'UA-36028965-1',
       "disqus_shortname" : 'labibi',
       "github_base_account" : 'ctb',
       "github_project" : 'labibi',
     }
