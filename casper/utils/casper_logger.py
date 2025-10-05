def clog(msg, title="CASPER LOGGING"):
    tlen = len(f"-------------------{title}-------------------")
    print(f'''
-------------------{title}-------------------
{msg}
{"-"*tlen}
''')
    
def rclog(msg, title="CASPER LOGGING"):
    tlen = len(f"-------------------{title}-------------------")
    return(f'''
-------------------{title}-------------------
{msg}
{"-"*tlen}
''')