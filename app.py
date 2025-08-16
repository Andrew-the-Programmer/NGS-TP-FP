import tkinter as tk
from collections import defaultdict
from tkinter import *
from tkinter import BOTH, HORIZONTAL, Frame, StringVar

from keys import Keys, format_percent
from tpfp import getTPFP


class Window(Frame):
    def __init__(self, master=None):
        Frame.__init__(self, master)
        self.master = master

        self.err_label = None

        self.params_coverage = (
            (Keys.VAF, "VAF (%)", "3"),
            (Keys.VR, "Variant reads (number)", ""),
            (Keys.COV, "Coverage depth (number)", "500"),
            (Keys.ERR, "Sequencing error (%)", "1"),
            (Keys.MVR, "Minimum variant reads (number)", "10"),
        )

        self.results_coverage = (
            (Keys.ERR_MSG, "Error message"),
            (
                Keys.FP,
                lambda v: f"Probability of false positive: {format_percent(v)}",
            ),
            # (
            #     Keys.TP,
            #     lambda v: f"Probability of true positive: {format_percent(v)}",
            # ),
        )

        self.entries: dict[str, tk.Entry] = {}
        self.entries_list: list[tk.Entry] = []
        self.results_strvars: defaultdict[str, StringVar] = defaultdict(StringVar)
        self.results_labels: dict[str, tk.Label] = {}

        self.initWindow()

    def entry_next(self, r):
        def callback(event):
            nexti = (r + 1) % len(self.params_coverage)
            entry = self.entries_list[nexti]
            entry.focus()
            entry.selection_adjust(tk.END)

        return callback

    def entry_prev(self, r):
        def callback(event):
            nexti = (r - 1) % len(self.params_coverage)
            entry = self.entries_list[nexti]
            entry.focus()
            entry.selection_adjust(tk.END)

        return callback

    def initWindow(self):
        self.master.title("OLGEN Coverage Limit")
        self.pack(fill=BOTH, expand=1)

        parent_panel = tk.PanedWindow(self, orient=HORIZONTAL)
        parent_panel.pack(fill=BOTH, expand=0)

        coverage_frame = tk.LabelFrame(parent_panel)
        coverage_frame.pack(fill=BOTH, expand=1)
        # parent_panel.add(coverage_frame)

        coverage_item_frame = Frame(coverage_frame, pady=10)
        coverage_item_frame.pack()

        r = 0
        for param in self.params_coverage:
            item_label = tk.Label(coverage_item_frame, text=param[1])
            item_label.grid(row=r, column=0)
            item_entry = tk.Entry(coverage_item_frame)
            item_entry.grid(row=r, column=1)
            item_entry.bind("<Return>", self.computeCoverage)
            item_entry.bind("<Down>", self.entry_next(r))
            item_entry.bind("<Up>", self.entry_prev(r))
            item_entry.insert(tk.END, param[2])
            if r == 0:
                item_entry.focus()

            self.entries[param[0]] = item_entry
            self.entries_list.append(item_entry)

            r += 1

        compute_button = tk.Button(
            coverage_frame, text="Compute errors", command=self.computeCoverage
        )
        compute_button.pack(side=tk.BOTTOM)

        for param in self.results_coverage[::-1]:
            result_label = tk.Label(
                coverage_frame,
                textvariable=self.results_strvars[param[0]],
                font=("Helvetica", 10, "bold"),
            )
            result_label.pack(side=tk.BOTTOM)
            self.results_labels[param[0]] = result_label

    def errMessage(self, msg):
        self.results_strvars[Keys.ERR_MSG].set(msg)
        self.results_labels[Keys.ERR_MSG].configure(fg="red")

    def warnMessage(self, msg):
        self.results_strvars[Keys.ERR_MSG].set(msg)
        self.results_labels[Keys.ERR_MSG].configure(fg="orange")

    def computeCoverage(self, component=None):
        self.results_labels[Keys.ERR_MSG].configure(fg="black")

        for result in self.results_strvars.values():
            result.set("")

        vaf = 0
        cov = 0

        vaf = self.entries[Keys.VAF].get()
        vr = self.entries[Keys.VR].get()
        cov = self.entries[Keys.COV].get()
        err = self.entries[Keys.ERR].get()
        mvr = self.entries[Keys.MVR].get()

        try:
            cov = int(cov)
        except ValueError:
            self.errMessage("Cannot convert coverage to int.")
            return

        if cov < 0:
            self.errMessage("Coverage must be non-negative.")
            return

        try:
            err = float(self.entries[Keys.ERR].get())
            err /= 100
        except ValueError:
            self.errMessage("Cannot convert sequencing error to float.")
            return

        if vaf == "" and vr == "":
            self.errMessage("VAF or variant reads must be specified.")
            return
        if vaf != "" and vr != "":
            self.errMessage(
                "VAF and variant reads cannot be specified at the same time."
            )
            return
        if vaf == "":
            try:
                vr = int(vr)
            except ValueError:
                self.errMessage("Cannot convert variant reads to int.")
                return
            vaf = vr / cov
        else:
            try:
                vaf = float(vaf)
                vaf /= 100
            except ValueError:
                self.errMessage("Cannot convert variant allele frequency to float.")
                return
            vr = round(cov * vaf)

        try:
            mvr = int(mvr)
        except ValueError:
            self.errMessage("Cannot convert minimum variant reads to int.")
            return

        try:
            results: dict[str] = getTPFP(cov, vr, err, mvr)
        except ValueError as e:
            self.errMessage(str(e))
            return

        prefixes = dict(self.results_coverage)

        for k, v in results.items():
            msg = prefixes[k]
            self.results_strvars[k].set(msg(v))


root = Tk()
app = Window(root)
root.mainloop()
