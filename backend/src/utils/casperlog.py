def clog(msg, title="CASPER LOGGING"):
    tlen = len(f"-------------------{title}-------------------")
    print(f'''
-------------------{title}-------------------
{msg}
{"-"*tlen}
''')