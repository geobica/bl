import time
	
def format_duration(seconds):
    seconds = max(0,int(seconds))
    m,s = divmod(seconds,60)
    return f'{m}:{s:02d}'

def progress_bar(prefix,frac,start_time,done=False):
    if done:
        print('\r'+' '*100+'\r',end='',flush=True)
        return
    elapsed = time.time()-start_time
    frac = min(max(frac,1e-6),1.0)
    remaining = elapsed*(1-frac)/frac
    msg = f'{prefix} [{int(frac*100)}% done, {format_duration(elapsed)} elapsed so far, estimated {format_duration(remaining)} remaining]'
    print('\r'+msg.ljust(90),end='',flush=True)