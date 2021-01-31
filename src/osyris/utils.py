


# def value_to_string(val, precision=3):
#     if (not isinstance(val, float)) or (val == 0):
#         text = str(val)
#     elif (abs(val) >= 10.0**(precision+1)) or \
#          (abs(val) <= 10.0**(-precision-1)):
#         text = "{val:.{prec}e}".format(val=val, prec=precision)
#     else:
#         text = "{}".format(val)
#         if len(text) > precision + 2 + (text[0] == '-'):
#             text = "{val:.{prec}f}".format(val=val, prec=precision)
#     return text
