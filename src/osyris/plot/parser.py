import numpy as np
from matplotlib.colors import LogNorm, Normalize


def get_norm(norm=None, vmin=None, vmax=None):
    if norm == "log":
        if vmin is not None:
            vmin = np.log10(vmin)
        if vmax is not None:
            vmax = np.log10(vmax)
        norm = LogNorm(vmin=vmin, vmax=vmax)
    else:
        norm = Normalize(vmin=vmin, vmax=vmax)
    return norm



def parse_layer(entry, mode=None, norm=None, vmin=None, vmax=None, operation=None, **kwargs):
    # mode = get_mode(mode)
    
    # print("NORM IS", norm)

    if isinstance(entry, dict):
        params = {key: entry[key] for key in set(entry.keys()) - set(["data", "mode", "operation"])}
        if "norm" not in params:
            params["norm"] = norm
        if "vmin" in params:
            vmin = params["vmin"]
            del params["vmin"]
        if "vmax" in params:
            vmax = params["vmax"]
            del params["vmax"]

        params["norm"] = get_norm(norm=params["norm"],
            vmin=vmin, vmax=vmax)

        print(vmin, vmax, params["norm"].vmin, params["norm"].vmax)

        for key, arg in kwargs.items():
            if key not in params:
                params[key] = arg
        # print(params)
        # params = {key: entry[key] for key in set(entry.keys() - set(["data"]))}

        settings = {}
        for key in ["mode", "operation"]:
            settings[key] = entry[key] if key in entry else eval(key)
            # if key in entry:
            #     settings[key] = entry["mode"]
            # else:
            #     layer_mode = mode
        # if settings["mode"] == "image":
        #     settings["mode"] = "imshow"
        # settings["mode"] = get_mode(settings["mode"])

        return entry["data"], settings, params
        # params = {key: entry[key] for key in set(entry.keys() - set(["data"]))}
        # return entry["data"], kwargs
        # ret
    else:
        params = {"norm": get_norm(norm=norm, vmin=vmin, vmax=vmax)}
        settings = {"mode": mode, "operation": operation}
        params.update(kwargs)
        # if settings["mode"] == "image":
        #     settings["mode"] = "imshow"
        return entry, settings, params
