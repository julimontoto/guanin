from . import guanin


def main():
    args = guanin.argParser()
    guanin.runQCview(args)
    guanin.runQCfilter(args)
    guanin.technorm(args)
    guanin.contnorm(args)
    guanin.evalnorm(args)
