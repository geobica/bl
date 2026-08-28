import curses

class MenuItem:
    def __init__(self,label,value,color=None,annotation=None):
        self.label = label
        self.value = value
        self.color = color
        self.annotation = annotation

class MenuSection:
    def __init__(self,header,items=None):
        self.header = header
        self.items = items

def _row_len(section):
    return len(section.items) if section.items is not None else 1

def _color_pair_for(color):
    return {
        None: 0,
        'light_gray': 0,
        'dark_gray': 4,
        'red': 2,
        'orange': 3,
    }.get(color,0)

def _attr_for(color,current,colors_ok):
    if not colors_ok:
        if current:
            return curses.A_REVERSE | curses.A_BOLD
        return curses.A_DIM if color=='dark_gray' else curses.A_NORMAL
    if current:
        return curses.color_pair(1) | curses.A_BOLD
    pair = _color_pair_for(color)
    attr = curses.color_pair(pair) if pair else curses.A_NORMAL
    if color=='dark_gray':
        attr |= curses.A_DIM
    return attr

def _setup_colors():
    try:
        curses.start_color()
        if getattr(curses,'COLORS',0)<=0:
            return False
    except Exception:
        return False
    try:
        curses.use_default_colors()
        bg = -1
    except Exception:
        bg = curses.COLOR_BLACK
    try:
        curses.init_pair(1,curses.COLOR_GREEN,bg)
        curses.init_pair(2,curses.COLOR_RED,bg)
        curses.init_pair(3,curses.COLOR_YELLOW,bg)
        curses.init_pair(4,curses.COLOR_WHITE,bg)
    except Exception:
        return False
    return True

def run_menu(sections,question,enable_overwrite_toggle=False):
    if not sections:
        return None,False
    row = next((i for i,section in enumerate(sections) if _row_len(section)>0),0)
    col = 0
    colors_ok = False
    overwrite = False

    def draw(stdscr):
        stdscr.erase()
        max_y,max_x = stdscr.getmaxyx()
        out_row = 0
        stdscr.addstr(out_row,0,question[:max_x-1]); out_row += 1
        current_label = None
        current_annotation = None
        for si,section in enumerate(sections):
            marker = '>' if si==row else ' '
            line = f'  {marker} {section.header}'
            try:
                stdscr.addstr(out_row,0,line[:max_x-1])
            except curses.error:
                pass
            out_row += 1
            if section.items is None:
                if si==row:
                    current_label = section.header
                    current_annotation = None
                continue
            prefix = '      '
            try:
                stdscr.addstr(out_row,0,prefix)
            except curses.error:
                pass
            c = len(prefix)
            item_sep = '   '
            for ii,item in enumerate(section.items):
                is_current = (si==row and ii==col)
                attr = _attr_for(item.color,is_current,colors_ok)
                text = item.label
                try:
                    stdscr.addstr(out_row,c,text[:max(0,max_x-1-c)],attr)
                except curses.error:
                    pass
                c += len(text)
                try:
                    stdscr.addstr(out_row,c,item_sep[:max(0,max_x-1-c)])
                except curses.error:
                    pass
                c += len(item_sep)
                if is_current:
                    current_label = item.label
                    current_annotation = item.annotation
            out_row += 1
        footer_row = out_row
        in_progress = f' ({current_annotation})' if current_annotation else ''
        overwrite_tag = ' [overwrite]' if overwrite else ''
        footer = f'Press [Enter] to select {current_label}{overwrite_tag}{in_progress}, use arrow keys to change selection.'
        try:
            stdscr.addstr(footer_row,0,footer[:max_x-1])
        except curses.error:
            pass
        stdscr.refresh()

    def loop(stdscr):
        nonlocal row,col,colors_ok,overwrite
        try:
            curses.curs_set(0)
        except Exception:
            pass
        colors_ok = _setup_colors()
        stdscr.keypad(True)
        while True:
            draw(stdscr)
            key = stdscr.getch()
            n = len(sections)
            if key in (curses.KEY_LEFT,ord('h')):
                width = _row_len(sections[row])
                if width>1:
                    col = (col-1)%width
            elif key in (curses.KEY_RIGHT,ord('l')):
                width = _row_len(sections[row])
                if width>1:
                    col = (col+1)%width
            elif key in (curses.KEY_UP,ord('k')):
                row = (row-1)%n
                col = min(col,_row_len(sections[row])-1)
            elif key in (curses.KEY_DOWN,ord('j')):
                row = (row+1)%n
                col = min(col,_row_len(sections[row])-1)
            elif enable_overwrite_toggle and key==ord('o'):
                overwrite = not overwrite
            elif key in (10,13,curses.KEY_ENTER):
                section = sections[row]
                if section.items is None:
                    return section.header,None,overwrite
                return section.header,section.items[col].value,overwrite
            elif key in (ord('q'),27):
                return None,None,False

    header,value,overwrite_result = curses.wrapper(loop)
    return value,overwrite_result
