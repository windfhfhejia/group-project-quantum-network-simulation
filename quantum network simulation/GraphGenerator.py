# -*- coding: utf-8 -*-
"""
Created on Wed May  8 09:47:36 2019

@author: xv18766
"""

import csv
import matplotlib.pyplot as plt
import os


#def miniplot_from_csv(
#        file_prefix, xvalue, yvalue)


def plot_from_csv(file_prefix, xvalue, yvalues, xlabel=None, ylabel=None,
                  title=None, suptitle=None, legend=None,
                  xscale=['linear'], yscale=['linear'], ypositive=False,
                  invert_xaxis=False, invert_yaxis=False,
                  output_file_prefix=None, formats=['svg'],
                  show_data_source=False, suppress_output=True):
    filename = file_prefix + '.csv'
    if legend == None:
        legend = {}
        for yvalue in yvalues:
            legend[yvalue] = yvalue
    with open(filename, 'r', newline = '') as csvfile:
        reader = csv.DictReader(csvfile)
        data = {xvalue: []}
        for yvalue in yvalues:
            data[yvalue] = []
        for row in reader:
            for key in data:
                data[key].append(float(row[key]))
        for yvalue in yvalues:
            plt.plot(xvalue, yvalue, data=data, label=legend[yvalue])
        plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
        if not xlabel == None:
            plt.xlabel(*xlabel)
        if not ylabel == None:
            plt.ylabel(*ylabel)
        if not title == None:
            plt.title(*title)
        if not suptitle == None:
            plt.suptitle(*suptitle)
        if show_data_source == True:
            plt.figtext(x=0, y=-0.15, s='Generated from ' + filename +
                        ' with data from columns:\n'
                        + ',\n'.join([xvalue, *yvalues]),
                        verticalalignment="top")
        plt.xscale(*xscale)
        plt.yscale(*yscale)
        if ypositive == True:
            plt.ylim(ymin=0)
        if invert_xaxis == True:
            plt.gca().invert_xaxis()
        if invert_yaxis == True:
            plt.gca().invert_yaxis()
        plt.grid(False)
        if output_file_prefix == 'default':
            for format in formats:
                plt.savefig(file_prefix + '.' + format,
                            format=format, bbox_inches='tight')
        elif output_file_prefix == None:
            pass
        else:
            for format in formats:
                plt.savefig(output_file_prefix + '.' + format,
                            format=format, bbox_inches='tight')
        plt.show()
        plt.close()
    if suppress_output == True:
        return
    else:
        return data


def plot_twinx_from_csv(
                file_prefix, xvalue, yvalues1, yvalues2,
                xlabel=None, ylabel1=None, ylabel2=None,
                title=None, suptitle=None, legend1=None, legend2=None,
                legends_on_both_sides=True,
                xscale=['linear'], yscale1=['linear'], yscale2=['linear'],
                ypositive1=False, ypositive2=False,
                invert_xaxis=False, invert_yaxis=False,
                output_file_prefix=None, formats=['svg'],
                show_data_source=False, suppress_output=True):
    filename = file_prefix + '.csv'
    if legend1 == None:
        legend1 = {}
        for yvalue in yvalues1:
            legend1[yvalue] = yvalue
    if legend2 == None:
        legend2 = {}
        for yvalue in yvalues2:
            legend1[yvalue] = yvalue
    with open(filename, 'r', newline = '') as csvfile:
        reader = csv.DictReader(csvfile)
        data = {xvalue: []}
        for yvalue in yvalues1:
            data[yvalue] = []
        for yvalue in yvalues2:
            data[yvalue] = []
        for row in reader:
            for key in data:
                data[key].append(float(row[key]))
        ax1 = plt.gca()
        ax2 = ax1.twinx()
        for yvalue in yvalues1:
            ax1.plot(xvalue, yvalue, data=data, label=legend1[yvalue])
        for yvalue in yvalues2:
            ax2.plot(xvalue, yvalue, '--', data=data, label=legend2[yvalue])
        if legends_on_both_sides == True:
            ax1.legend(bbox_to_anchor=(-0.14, 1), loc="upper right")
        else:
            ax1.legend(bbox_to_anchor=(1.14, 0), loc="lower left")
        ax2.legend(bbox_to_anchor=(1.14, 1), loc="upper left")
        if not xlabel == None:
            ax1.set_xlabel(*xlabel)
        if not ylabel1 == None:
            ax1.set_ylabel(*ylabel1)
        if not ylabel2 == None:
            ax2.set_ylabel(*ylabel2)
        if not title == None:
            plt.title(*title)
        if not suptitle == None:
            plt.suptitle(*suptitle)
        if show_data_source == True:
            plt.figtext(x=0, y=-0.15, s='Generated from ' + filename +
                        ' with data from columns:\n'
                        + ',\n'.join([xvalue, *yvalues1, *yvalues2]),
                        verticalalignment="top")
        ax1.set_xscale(*xscale)
        ax1.set_yscale(*yscale1)
        ax2.set_yscale(*yscale2)
        if ypositive1 == True:
            ax1.set_ylim(ymin=0)
        if ypositive2 == True:
            ax2.set_ylim(ymin=0)
        if invert_xaxis == True:
            plt.gca().invert_xaxis()
        if invert_yaxis == True:
            plt.gca().invert_yaxis()
        plt.grid(False)
        if output_file_prefix == 'default':
            for format in formats:
                plt.savefig(file_prefix + '.' + format,
                            format=format, bbox_inches='tight')
        elif output_file_prefix == None:
            pass
        else:
            for format in formats:
                plt.savefig(output_file_prefix + '.' + format,
                            format=format, bbox_inches='tight')
        if output_file_prefix == 'default':
            plt.savefig(file_prefix + '.' + format,
                        format=format, bbox_inches='tight')
        elif output_file_prefix == None:
            pass
        else:
            plt.savefig(output_file_prefix + '.' + format,
                        format=format, bbox_inches='tight')
        plt.show()
        plt.close()
    if suppress_output == True:
        return
    else:
        return data

def multiplot_from_csv(directory,
        file_prefixes, xvalue, yvalue,  linestyle='-',
        xlabel=None, ylabel=None,
        title=None, suptitle=None, legend=None,
        xscale=['linear'], yscale=['linear'], ypositive=False,
        invert_xaxis=False, invert_yaxis=False,
        output_file_prefix=None, formats=['svg'],
        show_data_source=False, suppress_output=True):
    os.chdir(directory)
    filenames = {file + '.csv': file for file in file_prefixes}
    if legend == None:
        legend = filenames
    for filename in filenames:
        with open(filename, 'r', newline = '') as csvfile:
            reader = csv.DictReader(csvfile)
            data = {xvalue: [], yvalue: []}
            for row in reader:
                for key in data:
                    data[key].append(float(row[key]))
            plt.plot(xvalue, yvalue, linestyle,
                     data=data, label=legend[filename])
    plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
    if not xlabel == None:
        plt.xlabel(*xlabel)
    if not ylabel == None:
        plt.ylabel(*ylabel)
    if not title == None:
        plt.title(*title)
    if not suptitle == None:
        plt.suptitle(*suptitle)
    if show_data_source == True:
        plt.figtext(x=0, y=-0.15, s='Generated from '
                    + ',\n'.join(filenames) +
                    ' with data from column:\n'
                    + ',\n'.join([xvalue, yvalue]),
                    verticalalignment="top")
    plt.xscale(*xscale)
    plt.yscale(*yscale)
    if ypositive == True:
        plt.ylim(ymin=0)
    if invert_xaxis == True:
        plt.gca().invert_xaxis()
    if invert_yaxis == True:
        plt.gca().invert_yaxis()
    plt.grid(False)
    if output_file_prefix == 'default':
        for format in formats:
            plt.savefig('multiplot ' + xvalue + ' vs ' + yvalue
                        + '.' + format,
                        format=format, bbox_inches='tight')
    elif output_file_prefix == None:
        pass
    else:
        for format in formats:
            plt.savefig(output_file_prefix + '.' + format,
                        format=format, bbox_inches='tight')
    plt.show()
    plt.close()
    if suppress_output == True:
        return
    else:
        return data