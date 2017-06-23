import numpy as np
import config_osiris as conf

class OsirisCommon:
     
    def __init__(self):
                
        return

    #=======================================================================================
    # The new field function is used to create a new data field. Say you want to take the
    # log of the density. You create a new field by calling:
    # mydata.new_field(name="log_rho",operation="np.log10(density)",unit="g/cm3",label="log(Density)")
    # The operation string is then evaluated using the 'eval' function.
    #=======================================================================================
    def new_field(self,name,operation="",unit="",label="",verbose=True,values=[]):
        
        if (len(operation) == 0) and (len(values) > 0):
            new_data = values
            op_parsed = operation
            depth = -1
        else:
            [op_parsed,depth] = self.parse_operation(operation)
            try:
                new_data = eval(op_parsed)
            except NameError:
                if verbose:
                    print("Error parsing operation when trying to create variable: "+name)
                    print("The attempted operation was: "+op_parsed)
                return
        self.data[name] = dict()
        self.data[name]["values"   ] = new_data
        self.data[name]["unit"     ] = unit
        self.data[name]["label"    ] = label
        self.data[name]["operation"] = op_parsed
        self.data[name]["depth"    ] = depth+1
        
        return
    
    #=======================================================================================
    # The operation parser converts an operation string into an expression which contains
    # variables from the data dictionary. If a name from the variable list, e.g. "density",
    # is found in the operation, it is replaced by self.data["density"]["values"] so that it
    # can be properly evaluated by the 'eval' function in the 'new_field' function.
    #=======================================================================================
    def parse_operation(self,operation):
        
        max_depth = 0
        # Add space before and after to make it easier when searching for characters before
        # and after
        expression = " "+operation+" "
        # Sort the list of variable keys in the order of the longest to the shortest.
        # This guards against replacing 'B' inside 'logB' for example.
        key_list = sorted(self.data.keys(),key=lambda x:len(x),reverse=True)
        # For replacing, we need to create a list of hash keys to replace on instance at a
        # time
        hashkeys  = dict()
        hashcount = 0
        for key in key_list:
            # Search for all instances in string
            loop = True
            loc = 0
            while loop:
                loc = expression.find(key,loc)
                if loc == -1:
                    loop = False
                else:
                    # Check character before and after. If they are either a letter or a '_'
                    # then the instance is actually part of another variable or function
                    # name.
                    char_before = expression[loc-1]
                    char_after  = expression[loc+len(key)]
                    bad_before = (char_before.isalpha() or (char_before == "_"))
                    bad_after = (char_after.isalpha() or (char_after == "_"))
                    hashcount += 1
                    if (not bad_before) and (not bad_after):
                        theHash = "#"+str(hashcount).zfill(5)+"#"
                        # Store the data key in the hash table
                        hashkeys[theHash] = "self.data[\""+key+"\"][\"values\"]"
                        expression = expression.replace(key,theHash,1)
                        max_depth = max(max_depth,self.data[key]["depth"])
                    else:
                        # Replace anyway to prevent from replacing "x" in "max("
                        theHash = "#"+str(hashcount).zfill(5)+"#"
                        hashkeys[theHash] = key
                        expression = expression.replace(key,theHash,1)
                    loc += 1
        # Now go through all the hashes in the table and build the final expression
        for theHash in hashkeys.keys():
            expression = expression.replace(theHash,hashkeys[theHash])
    
        return [expression,max_depth]
    
    #=======================================================================================
    # The function get returns the values of the selected variable
    #=======================================================================================
    def get(self,var):
        return self.data[var]["values"]
