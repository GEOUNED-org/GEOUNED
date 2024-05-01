def comments_write(name, MetaList):
    """Function to write in an independent file the comment strings"""

    outfile = open(name + "_comments.txt", "w", encoding="utf-8")
    for m in MetaList:
        outfile.write(m.Comments + "\n")

    outfile.close()


def summary_write(name, MetaList):

    outfile = open(name + "_Summary.txt", "w", encoding="utf-8")
    header = f"  Cell Id{'':5s}Mat Id{'':6s}Density{'':7s}Volume{'':5s}Comments\n"
    outfile.write(header)

    for m in MetaList:
        if m.Void or m.__id__ is None:
            continue
        index = m.__id__
        Vol = m.Volume * 1e-3
        if m.Material == 0:
            line = " {:>8d}{:3s}{:>8d}{:3s}{:11s}{:3s}{:11.4e}{:3s}{}\n".format(
                index, "", 0, "", "", "", Vol, "", m.Comments
            )
        else:
            if abs(m.Density) < 1e-2:
                line = " {:>8d}{:3s}{:>8d}{:3s}{:11.4e}{:3s}{:11.4e}{:3s}{}\n".format(
                    index, "", m.Material, "", m.Density, "", Vol, "", m.Comments
                )
            else:
                line = " {:>8d}{:3s}{:>8d}{:3s}{:11.7f}{:3s}{:11.4e}{:3s}{}\n".format(
                    index, "", m.Material, "", m.Density, "", Vol, "", m.Comments
                )

        outfile.write(line)
    outfile.close()
