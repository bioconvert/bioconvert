###########################################################################
# Bioconvert is a project to facilitate the interconversion               #
# of life science data from one format to another.                        #
#                                                                         #
# Copyright © 2018-2022  Institut Pasteur, Paris and CNRS.                #
#                                                                         #
# bioconvert is free software: you can redistribute it and/or modify      #
# it under the terms of the GNU General Public License as published by    #
# the Free Software Foundation, either version 3 of the License, or       #
# (at your option) any later version.                                     #
#                                                                         #
# bioconvert is distributed in the hope that it will be useful,           #
# but WITHOUT ANY WARRANTY; without even the implied warranty of          #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
# GNU General Public License for more details.                            #
#                                                                         #
# You should have received a copy of the GNU General Public License       #
# along with this program (COPYING file).                                 #
# If not, see <http://www.gnu.org/licenses/>.                             #
#                                                                         #
# Repository: https://github.com/bioconvert/bioconvert                    #
# Documentation: http://bioconvert.readthedocs.io                         #
###########################################################################
"""Network tools to manipulate the graph of conversion"""
from os import environ

from bioconvert.core.registry import Registry

import colorlog

_log = colorlog.getLogger(__name__)


__all__ = ["create_graph", "get_conversions_wrapped", "create_graph_for_cytoscape"]


def create_graph(
    filename, layout="dot", use_singularity=False, color_for_disabled_converter="red", include_subgraph=False
):
    """

    :param filename: should end in .png or .svg or .dot

    If extension is .dot, only the dot file is created without annotations.
    This is useful if you have issues installing graphviz.
    If so, under Linux you could use our singularity container
    see github.com/cokelaer/graphviz4all

    """

    rr = Registry()
    try:
        if filename.endswith(".dot") or use_singularity is True:
            raise Exception()
        # local import because optional for bioconvert
        _log.info("Switching to pygraphviz")
        from pygraphviz import AGraph

        dg = AGraph(directed=True)

        url = "https://bioconvert.readthedocs.io/en/main/formats.html#{}"

        for a, b, s in rr.get_all_conversions():
            if len(a) == 1 and len(b) == 1:

                dg.add_node(
                    a[0],
                    shape="rectangle",
                    style="filled",
                    url=url.format(a[0].upper()),
                )
                dg.add_node(
                    b[0],
                    shape="rectangle",
                    style="filled",
                    url=url.format(b[0].upper()),
                )
                dg.add_edge(a[0], b[0], alpha=0.5, color="black" if s else color_for_disabled_converter, minlen=1)
            else:
                and_node = "_".join(a) + "_and_" + "_".join(b)

                dg.add_node(
                    and_node,
                    label="",
                    fillcolor="black",
                    width=0.1,
                    height=0.1,
                    styled="filled",
                    fixedsize=True,
                    shape="circle",
                )

                for this in a:
                    dg.add_edge(
                        this,
                        and_node,
                        color="black" if s else color_for_disabled_converter,
                    )

                for this in b:
                    dg.add_edge(
                        and_node,
                        this,
                        color="black" if s else color_for_disabled_converter,
                    )

        for name in dg.nodes():
            if dg.degree(name) < 5:
                dg.get_node(name).attr["fillcolor"] = "white"
            elif dg.degree(name) < 10:
                # yellow
                dg.get_node(name).attr["fillcolor"] = "yellow"
            elif dg.degree(name) < 20:
                # orange
                dg.get_node(name).attr["fillcolor"] = "orange"
            else:
                # red
                dg.get_node(name).attr["fillcolor"] = "red"

        if include_subgraph:
            # sequencing
            with dg.subgraph(name="cluster_sequencing", shape="circle") as c:
                c.graph_attr.update(style="filled", color="lightgrey", fillcolor="#A569BD", label="Sequencing")
                c.add_nodes_from(
                    [
                        "FASTQ",
                        "FASTA",
                        "SRA",
                        "FAA",
                        "AGP",
                        "TWOBIT",
                        "ABI",
                        "QUAL",
                        "FASTA_QUAL_and_FASTQ",
                        "FASTQ_and_FASTA_QUAL",
                        "FASTA_and_FASTA_AGP",
                    ]
                )

            # alignment
            with dg.subgraph(name="cluster_alignment") as c:
                c.graph_attr.update(
                    style="filled", color="lightgrey", shape="circle", fillcolor="#D2B4DE", label="Alignment"
                )
                c.add_nodes_from(["SAM", "BAM", "CRAM", "PAF", "MAF"])

            # phylogney
            with dg.subgraph(name="cluster_phylo") as c:
                c.graph_attr.update(
                    style="filled", color="lightgrey", shape="box", fillcolor="#BB8FCE", label="Phylogney"
                )
                c.add_nodes_from(["NEXUS", "PHYLOXML", "CLUSTAL", "NEWICK", "PHYLIP", "STOCKHOLM"])

            # variant
            with dg.subgraph(name="cluster_variant") as c:
                c.graph_attr.update(
                    style="filled", color="lightgrey", shape="box", fillcolor="#F4ECF7", label="Variant"
                )
                c.add_nodes_from(["VCF", "BCF", "PLINK", "BPLINK"])

            # annotation
            with dg.subgraph(name="cluster_annotation") as c:
                c.graph_attr.update(
                    style="filled", color="lightgrey", shape="box", fillcolor="#D2B4DE", label="Annotation"
                )
                c.add_nodes_from(["GENBANK", "EMBL", "GFF3", "GFF2"])

            # compression
            with dg.subgraph(name="cluster_comp") as c:
                c.graph_attr.update(
                    style="filled", color="lightgrey", shape="box", fillcolor="#E9DAEF", label="Compression"
                )
                c.add_nodes_from(["GZ", "DSRC", "BZ2"])

            # coverage
            with dg.subgraph(name="cluster_cov") as c:

                c.graph_attr.update(
                    style="filled", color="lightgrey", shape="box", fillcolor="#E9DAEF", label="Coverage"
                )
                c.add_nodes_from(["COV", "BED", "BEDGRAPH", "WIG", "BIGWIG", "WIGGLE", "BIGBED"])

            # assembly
            with dg.subgraph(name="cluster_ass") as c:

                c.graph_attr.update(
                    style="filled", color="lightgrey", shape="box", fillcolor="#BB8FCE", label="Assembly"
                )
                c.add_nodes_from(["GFA", "SCF"])

        dg.layout(layout)
        dg.draw(filename)
        dg.write("conversion.dot")
        return dg

    except Exception as e:
        print(e)
        _log.info("Switching to local dot and singularity")
        dot = """
strict digraph{
    node [label="\\N"];

    """
        nodes = set([item for items in rr.get_all_conversions() for item in items[0:1][0]])
        for node in nodes:
            dot += '"{}";\n'.format(node)
        for a, b, s in rr.get_all_conversions():
            dot += '"{}" -> "{}";\n'.format(a[0], b[0])
        dot += "}\n"

        from bioconvert import shell
        from easydev import TempFile

        dotfile = TempFile(suffix=".dot")
        with open(dotfile.name, "w") as fout:
            fout.write(dot)

        dotpath = ""
        if use_singularity:
            from bioconvert.core.downloader import download_singularity_image

            try:
                singfile = download_singularity_image(
                    "graphviz.simg",
                    "shub://cokelaer/graphviz4all:v1",
                    "4288088d91c848e5e3a327282a1ab3d1",
                )

            except Exception:
                print(
                    "Warning ! Singularity must be installed if you want to you used it ! Switching to local graphviz executable if available"
                )

            else:
                dotpath = "singularity run {} ".format(singfile)
                on_rtd = environ.get("READTHEDOCS", None) == "True"
                if on_rtd:
                    dotpath = ""

        ext = filename.rsplit(".", 1)[1]
        cmd = "{}dot -T{} {} -o {}".format(dotpath, ext, dotfile.name, filename)
        try:
            shell(cmd)
        except Exception:
            import os

            os.system(cmd)


def get_conversions_wrapped(registry, all_conversions=False):
    if all_conversions:
        for i, o, s in registry.get_all_conversions():
            yield i, o, s
    else:
        for i, o in registry.get_conversions():
            yield i, o, True


def create_graph_for_cytoscape(all_converter=False):
    """

    :param all_converter: use all converters or only the ones
        available in the current installation
    :return:
    """
    from bioconvert.core.registry import Registry

    registry = Registry()
    graph_nodes = []
    graph_edges = []
    graph = {
        "data": {
            "selected": False,
        },
        "elements": {
            "nodes": graph_nodes,
            "edges": graph_edges,
        },
    }
    nodes = {}

    def get_or_create(fmt):
        try:
            return nodes[fmt]
        except:
            ret = {
                "data": {
                    "id": "n" + str(len(nodes)),
                    "name": fmt,
                }
            }
            nodes[fmt] = ret
            graph_nodes.append(ret)
            return ret

    for i, o, _ in get_conversions_wrapped(registry, all_converter):
        i_as_node = get_or_create(i)
        o_as_node = get_or_create(o)
        graph_edges.append(
            {
                "data": {
                    "id": "e" + str(len(graph_edges)),
                    "source": i_as_node["data"]["id"],
                    "target": o_as_node["data"]["id"],
                }
            }
        )

    return graph
