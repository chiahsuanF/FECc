#!/usr/bin/env python3

import os
import logging

import matplotlib

# Required when running on a cluster without a display.
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def strandless_node(node):
    """
    Remove the terminal strand symbol.

    Example:
        chr2:100-200+ -> chr2:100-200
        chr2:100-200- -> chr2:100-200
    """
    node = node.strip()

    if node.endswith("+") or node.endswith("-"):
        return node[:-1]

    return node


def node_strand(node):
    """
    Return the terminal strand symbol.
    """
    if node.endswith("+") or node.endswith("-"):
        return node[-1]

    return "."


def parse_fragment_length(node):
    """
    Calculate fragment length from a coordinate node.

    Example:
        chr2:100-200+ -> 101 bp
    """
    region = strandless_node(node)

    chromosome, coordinates = region.rsplit(":", 1)
    start_text, end_text = coordinates.split("-", 1)

    start = int(start_text)
    end = int(end_text)

    return end - start + 1


def parse_bigcycle_path(path_text):
    """
    Parse a path written as:

        A -> B -> C -> A

    Semicolon-separated paths are also supported.
    """
    path_text = path_text.strip()

    if "->" in path_text:
        fields = path_text.split("->")
    else:
        fields = path_text.split(";")

    return [
        field.strip()
        for field in fields
        if field.strip()
    ]


def make_color_map(regions):
    """
    Assign one color to every unique strandless genomic region.
    """
    region_count = len(regions)

    if region_count <= 10:
        color_values = plt.cm.tab10(
            np.linspace(
                0,
                1,
                max(region_count, 1)
            )
        )

    elif region_count <= 20:
        color_values = plt.cm.tab20(
            np.linspace(
                0,
                1,
                region_count
            )
        )

    else:
        color_values = plt.cm.hsv(
            np.linspace(
                0,
                1,
                region_count,
                endpoint=False
            )
        )

    return {
        region: color_values[index]
        for index, region in enumerate(regions)
    }


def scale_marker_sizes(
    fragment_lengths,
    minimum_size=100,
    maximum_size=900
):
    """
    Scale marker area according to fragment length.

    Square-root scaling prevents very long fragments from completely
    dominating the figure.
    """
    lengths = np.asarray(
        fragment_lengths,
        dtype=float
    )

    if len(lengths) == 0:
        return np.asarray([])

    minimum_length = lengths.min()
    maximum_length = lengths.max()

    if minimum_length == maximum_length:
        return np.full(
            len(lengths),
            (minimum_size + maximum_size) / 2.0
        )

    normalized = (
        (lengths - minimum_length)
        /
        (maximum_length - minimum_length)
    )

    return (
        minimum_size
        + (maximum_size - minimum_size)
        * np.sqrt(normalized)
    )
    
def format_fragment_size(length):
    """
    Format a fragment length for display.

    Examples:
        157 bp
        135.8 kb
        1.42 Mb
    """
    if length >= 1000000:
        return "{:.2f} Mb".format(
            length / 1000000.0
        )

    if length >= 1000:
        return "{:.1f} kb".format(
            length / 1000.0
        )

    return "{} bp".format(length)


def plot_bigcycle_path(
    path_text,
    output_file,
    bigcycle_id=None,
    dpi=250,
    minimum_arc_degrees=5.0
):
    """
    Plot one eccDNA as a segmented circular ring.

    Each path occurrence is drawn as one colored arc.

    Arc color:
        Genomic fragment identity.

    Arc length:
        Proportional to the square root of fragment length. The square
        root transformation keeps very short fragments visible when
        very long fragments are also present.

    Arrowhead:
        Traversal direction.

    Labels:
        Chromosome, coordinates, strand, and fragment length.
    """
    nodes = parse_bigcycle_path(path_text)

    if len(nodes) < 2:
        logging.warning(
            "Cannot plot big cycle %s: fewer than two nodes",
            bigcycle_id
        )
        return None

    # The final node may repeat the first node only to close the cycle.
    is_closed = nodes[0] == nodes[-1]

    if is_closed:
        plot_nodes = nodes[:-1]
    else:
        plot_nodes = nodes

    if len(plot_nodes) < 2:
        logging.warning(
            "Cannot plot big cycle %s: fewer than two plotted fragments",
            bigcycle_id
        )
        return None

    unique_regions = list(
        dict.fromkeys(
            strandless_node(node)
            for node in plot_nodes
        )
    )

    color_map = make_color_map(unique_regions)

    fragment_lengths = np.asarray(
        [
            parse_fragment_length(node)
            for node in plot_nodes
        ],
        dtype=float
    )

    # Square-root scaling gives longer fragments larger arcs while
    # keeping small fragments visible.
    arc_weights = np.sqrt(fragment_lengths)

    node_count = len(plot_nodes)

    # Leave a small gap between neighboring fragments.
    gap_degrees = min(
        2.0,
        20.0 / float(max(node_count, 1))
    )

    available_degrees = (
        360.0
        - gap_degrees * node_count
    )

    if available_degrees <= 0:
        gap_degrees = 0.0
        available_degrees = 360.0

    raw_arc_degrees = (
        arc_weights
        / arc_weights.sum()
        * available_degrees
    )

    # Make tiny fragments visible.
    arc_degrees = np.maximum(
        raw_arc_degrees,
        minimum_arc_degrees
    )

    # Rescale after applying the minimum.
    arc_degrees = (
        arc_degrees
        / arc_degrees.sum()
        * available_degrees
    )

    figure_size = max(
        9,
        min(
            20,
            8 + node_count * 0.18
        )
    )

    figure, axis = plt.subplots(
        figsize=(figure_size, figure_size)
    )

    radius = 1.0
    line_width = 8.0

    # Begin at the top and move clockwise.
    current_angle = 90.0

    label_records = []

    for index, node in enumerate(plot_nodes):
        region = strandless_node(node)
        strand = node_strand(node)
        fragment_length = int(fragment_lengths[index])

        fragment_angle = float(
            arc_degrees[index]
        )

        start_angle = current_angle
        end_angle = (
            current_angle
            - fragment_angle
        )

        # Leave room for an arrowhead at the end of the arc.
        arrow_degrees = min(
            4.0,
            max(1.2, fragment_angle * 0.12)
        )

        arc_end_angle = (
            end_angle
            + arrow_degrees
        )

        # Draw the main colored arc.
        angle_values = np.linspace(
            np.radians(start_angle),
            np.radians(arc_end_angle),
            120
        )

        x_values = radius * np.cos(
            angle_values
        )

        y_values = radius * np.sin(
            angle_values
        )

        axis.plot(
            x_values,
            y_values,
            color=color_map[region],
            linewidth=line_width,
            solid_capstyle="round",
            zorder=2
        )

        # Draw an arrowhead following the direction of the circle.
        arrow_start_angle = np.radians(
            arc_end_angle
        )

        arrow_end_angle = np.radians(
            end_angle
        )

        arrow_start = (
            radius * np.cos(arrow_start_angle),
            radius * np.sin(arrow_start_angle)
        )

        arrow_end = (
            radius * np.cos(arrow_end_angle),
            radius * np.sin(arrow_end_angle)
        )

        axis.annotate(
            "",
            xy=arrow_end,
            xytext=arrow_start,
            arrowprops={
                "arrowstyle": "-|>",
                "color": color_map[region],
                "linewidth": 2.0,
                "mutation_scale": 16,
                "shrinkA": 0,
                "shrinkB": 0
            },
            zorder=3
        )

        # Label at the midpoint of the fragment arc.
        midpoint_angle = np.radians(
            (
                start_angle
                + end_angle
            ) / 2.0
        )

        label_radius = 1.27

        label_x = (
            label_radius
            * np.cos(midpoint_angle)
        )

        label_y = (
            label_radius
            * np.sin(midpoint_angle)
        )

        horizontal_alignment = (
            "left"
            if label_x >= 0
            else "right"
        )

        label_text = (
            "{}\n{}\nstrand {}\n{}".format(
                region.split(":", 1)[0],
                region.split(":", 1)[1],
                strand,
                format_fragment_size(fragment_length)
            )
        )

        label_records.append(
            (
                label_x,
                label_y,
                horizontal_alignment,
                label_text
            )
        )

        current_angle = end_angle - gap_degrees

    # Draw labels after drawing all arcs.
    for (
        label_x,
        label_y,
        horizontal_alignment,
        label_text
    ) in label_records:

        axis.text(
            label_x,
            label_y,
            label_text,
            horizontalalignment=horizontal_alignment,
            verticalalignment="center",
            fontsize=8
        )

    # Add central eccDNA label.
    center_lines = []

    if bigcycle_id is not None:
        center_lines.append(
            "FEC_eccDNA {}".format(bigcycle_id)
        )
    else:
        center_lines.append("eccDNA")

    center_lines.append(
        "{} fragments".format(node_count)
    )

    axis.text(
        0,
        0,
        "\n".join(center_lines),
        horizontalalignment="center",
        verticalalignment="center",
        fontsize=14,
        fontweight="bold"
    )

    # Add a color legend below the circular plot.
    for region in unique_regions:
        axis.plot(
            [],
            [],
            color=color_map[region],
            linewidth=6,
            label=region
        )

    axis.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.04),
        fontsize=7,
        frameon=False,
        ncol=min(3, len(unique_regions))
    )

    axis.set_xlim(-1.75, 1.75)
    axis.set_ylim(-1.65, 1.65)
    axis.set_aspect("equal")
    axis.axis("off")

    title = "Reconstructed eccDNA"
    """
    if source_cycles:
        title = "{}\nSource cycles: {}".format(
            title,
            source_cycles
        )
    """
    axis.set_title(
        title,
        fontsize=12,
        pad=18
    )

    figure.tight_layout()

    figure.savefig(
        output_file,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white"
    )

    plt.close(figure)

    logging.info(
        "Wrote circular eccDNA plot -> %s",
        output_file
    )

    return output_file


def plot_bigcycles_from_summary(
    summary_bigcycle_file,
    output_dir,
    dpi=200
):
    """
    Read summary_bigcycle.txt and create one PNG for each big cycle.
    """
    os.makedirs(
        output_dir,
        exist_ok=True
    )

    output_files = []

    with open(
        summary_bigcycle_file,
        "r"
    ) as handle:

        header_line = handle.readline().rstrip("\n")

        if not header_line:
            logging.warning(
                "Big-cycle summary is empty: %s",
                summary_bigcycle_file
            )
            return output_files

        header = header_line.split("\t")

        column_index = {
            name: index
            for index, name in enumerate(header)
        }

        required_columns = {
            "bigcycle_id",
            "path"
        }

        missing_columns = (
            required_columns
            - set(column_index)
        )

        if missing_columns:
            raise ValueError(
                "Missing required columns in {}: {}".format(
                    summary_bigcycle_file,
                    ",".join(sorted(missing_columns))
                )
            )

        for line_number, line in enumerate(
            handle,
            start=2
        ):
            line = line.rstrip("\n")

            if not line:
                continue

            fields = line.split("\t")

            if len(fields) < len(header):
                logging.warning(
                    "Skipping malformed line %d in %s: "
                    "expected %d fields, found %d",
                    line_number,
                    summary_bigcycle_file,
                    len(header),
                    len(fields)
                )
                continue

            bigcycle_id = fields[
                column_index["bigcycle_id"]
            ]

            path_text = fields[
                column_index["path"]
            ]

            output_file = os.path.join(
                output_dir,
                "bigcycle_{}.png".format(
                    bigcycle_id
                )
            )

            result = plot_bigcycle_path(
                path_text=path_text,
                output_file=output_file,
                bigcycle_id=bigcycle_id,
                dpi=dpi
            )

            if result is not None:
                output_files.append(result)

    logging.info(
        "Created %d big-cycle plots in %s",
        len(output_files),
        output_dir
    )

    return output_files
 
