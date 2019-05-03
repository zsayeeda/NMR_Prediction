def readChemicalShifts(file, c_pos, c_shift):
        #reader = None
        try:
                with open(file, 'r') as fp:
                    text = fp.readlines()
                parse_shifts = False
                for i in range(len(text)):
                    items = text[i].split("\t")
                    if parse_shifts:
                        if text[i] == "":
                            continue 
                        if "Table" in text[i]:
                            parse_shifts = False
                            continue 
                        print("The file with file path: {}".format(file))
                        if items[0].strip() == "No." | "M" in items[2]:
                            continue 
                        c_pos.append(items[1].strip())
                        if items[2] != None & items[2] != "":
                            c_shift.append(items[2].stip())
                        else:
                            c_shift.append(items[3].stip())
                        print("The file with file path: {}".format(file))
                    if len(items) == 0 | len(items) == 1:
                        continue 
                    if "Atom" in items[1].strip():
                        parse_shifts = True
        except FileNotFoundError as e:
            print(e)
        except IOError as e:
            print(e)
        finally:
            try:
                if fp != None:
                    fp.close()
            except IOError as e:
                print(e)