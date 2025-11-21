import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, Input, Output, State, callback_context, dash_table
import dash_bootstrap_components as dbc
from pathlib import Path
import sys
import json
import csv
import io
import base64
from typing import List, Optional, Dict, Any
from tqdm import tqdm
from colorama import init, Fore, Style
from datetime import datetime
import networkx as nx
from collections import defaultdict

init()


def get_colored_time():
    current_time = datetime.now().strftime("%H:%M:%S")
    return f"{Fore.CYAN}[{current_time}]{Style.RESET_ALL}"


class EnhancedNetworkVisualizer:
    def __init__(self, data_path: str, har_tf_path: Optional[str] = None):
        self.data_path = data_path
        self.har_tf_path = har_tf_path
        self.network_data = None
        self.har_tf_data = None
        self.subnetwork_data = None
        self._cache = {}  # 添加缓存机制
        self.load_data()

    def load_data(self):
        try:
            print(f"{get_colored_time()} Loading network data from {self.data_path}...")
            
            # 优化数据加载，只读取需要的列
            required_columns = ["TF", "Target", "Weight", "CellType", "Stage", "Region"]
            self.network_data = pd.read_csv(self.data_path, usecols=required_columns)

            print(f"{get_colored_time()} Cleaning and optimizing data...")
            
            # 批量处理数据清理，提高性能
            self.network_data = self.network_data.fillna({
                "Weight": 0,
                "TF": "Unknown_TF",
                "Target": "Unknown_Target", 
                "Stage": "Unknown_Stage",
                "Region": "Unknown_Region",
                "CellType": "Unknown_CellType"
            })
            
            # 过滤无效数据
            self.network_data = self.network_data[
                (self.network_data["Weight"] > 0) &
                (self.network_data["TF"] != "Unknown_TF") &
                (self.network_data["Target"] != "Unknown_Target")
            ].copy()
            
            # 优化数据类型
            self.network_data["Weight"] = self.network_data["Weight"].astype("float32")
            
            print(f"{get_colored_time()} Network data loaded successfully. Shape: {self.network_data.shape}")

            if self.har_tf_path:
                print(f"{get_colored_time()} Loading HAR-TF data from {self.har_tf_path}...")
                self.har_tf_data = pd.read_csv(self.har_tf_path)
                print(f"{get_colored_time()} HAR-TF data loaded successfully. Shape: {self.har_tf_data.shape}")
            else:
                print(f"{get_colored_time()} No HAR-TF data provided.")

        except Exception as e:
            print(f"{get_colored_time()} {Fore.RED}Error loading data: {str(e)}{Style.RESET_ALL}")
            sys.exit(1)

    def get_available_options(self):
        # 数据已经在load_data中过滤过了，直接使用
        return {
            "tfs": sorted(self.network_data["TF"].unique().tolist()),
            "targets": sorted(self.network_data["Target"].unique().tolist()),
            "stages": sorted(self.network_data["Stage"].unique().tolist()),
            "regions": sorted(self.network_data["Region"].unique().tolist()),
            "celltypes": sorted(self.network_data["CellType"].unique().tolist()),
            "har_tfs": sorted(self.har_tf_data["TF"].unique().tolist())
            if self.har_tf_data is not None
            else [],
        }

    def get_top_nodes(self, node_type: str, top_n: int, filters: Dict[str, Any] = None):
        # 数据已经在load_data中过滤过了，直接使用
        data = self.network_data.copy()

        if filters:
            for key, value in filters.items():
                if value and key in data.columns:
                    if isinstance(value, list):
                        data = data[data[key].isin(value)]
                    else:
                        data = data[data[key] == value]

        if node_type == "tf":
            top_nodes = (
                data.groupby("TF")["Weight"].sum().nlargest(top_n).index.tolist()
            )
        elif node_type == "target":
            top_nodes = (
                data.groupby("Target")["Weight"].sum().nlargest(top_n).index.tolist()
            )
        else:
            top_nodes = []

        return top_nodes

    def create_sankey_diagram(
        self,
        show_har: bool = True,
        top_n_tfs: int = 10,
        top_n_targets: int = 20,
        specific_tfs: List[str] = None,
        specific_targets: List[str] = None,
        specific_stages: List[str] = None,
        specific_regions: List[str] = None,
        specific_celltypes: List[str] = None,
        height: int = 600,
        width: int = 1000,
    ):
        # 数据已经在load_data中过滤过了，直接使用
        flow_data = self.network_data.copy()

        if specific_celltypes:
            flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]
        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]
        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]

        if len(flow_data) == 0:
            return None, "No data found for the specified filters"

        if specific_tfs:
            available_tfs = set(flow_data["TF"].unique())
            top_tfs = [tf for tf in specific_tfs if tf in available_tfs]
        else:
            top_tfs = self.get_top_nodes(
                "tf",
                top_n_tfs,
                {
                    "CellType": specific_celltypes,
                    "Region": specific_regions,
                    "Stage": specific_stages,
                },
            )

        if specific_targets:
            available_targets = set(flow_data["Target"].unique())
            top_targets = [
                target for target in specific_targets if target in available_targets
            ]
        else:
            flow_data_filtered = flow_data[flow_data["TF"].isin(top_tfs)]
            top_targets = self.get_top_nodes(
                "target",
                top_n_targets,
                {
                    "CellType": specific_celltypes,
                    "Region": specific_regions,
                    "Stage": specific_stages,
                },
            )

        flow_data = flow_data[
            flow_data["TF"].isin(top_tfs) & flow_data["Target"].isin(top_targets)
        ].copy()

        if len(flow_data) == 0:
            return None, "No connections found for selected nodes"

        levels_to_show = ["TF"]
        level_nodes = {
            "HAR": [],
            "TF": top_tfs,
            "Stage": sorted(flow_data["Stage"].unique()),
            "Region": sorted(flow_data["Region"].unique()),
            "CellType": sorted(flow_data["CellType"].unique()),
            "Target": top_targets,
        }

        har_to_tf_map = {}

        if show_har and self.har_tf_data is not None:
            levels_to_show.insert(0, "HAR")
            har_tf_groups = (
                self.har_tf_data[self.har_tf_data["TF"].isin(top_tfs)]
                .groupby("TF")["HAR"]
                .agg(lambda x: sorted(x))
                .reset_index()
            )

            har_nodes = []
            for _, row in har_tf_groups.iterrows():
                tf = row["TF"]
                har_list = row["HAR"]
                first_har = har_list[0]
                har_count = len(har_list)
                har_display = f"HARs: {first_har}...({har_count})"
                har_nodes.append(har_display)
                har_to_tf_map[har_display] = tf

            level_nodes["HAR"] = har_nodes

        if not any([specific_celltypes, specific_stages, specific_regions]):
            levels_to_show.extend(["Stage", "Region", "CellType"])
        else:
            if not specific_stages:
                levels_to_show.append("Stage")
            if not specific_regions:
                levels_to_show.append("Region")
            if not specific_celltypes:
                levels_to_show.append("CellType")

        levels_to_show.append("Target")

        nodes = []
        node_to_idx = {}
        current_idx = 0

        for level in levels_to_show:
            for node in level_nodes[level]:
                nodes.append(node)
                node_to_idx[node] = current_idx
                current_idx += 1

        sources = []
        targets = []
        values = []
        link_colors = []
        all_connections = []

        for i in range(len(levels_to_show) - 1):
            current_level = levels_to_show[i]
            next_level = levels_to_show[i + 1]

            if current_level == "HAR" and next_level == "TF":
                for har_group in level_nodes["HAR"]:
                    tf = har_to_tf_map[har_group]
                    if har_group in node_to_idx and tf in node_to_idx:
                        sources.append(node_to_idx[har_group])
                        targets.append(node_to_idx[tf])
                        tf_weight = flow_data[flow_data["TF"] == tf]["Weight"].sum()
                        values.append(float(tf_weight))
                        all_connections.append(tf)
            else:
                grouped_data = (
                    flow_data.groupby([current_level, next_level])["Weight"]
                    .sum()
                    .reset_index()
                )

                for _, row in grouped_data.iterrows():
                    if (
                        row[current_level] in node_to_idx
                        and row[next_level] in node_to_idx
                    ):
                        sources.append(node_to_idx[row[current_level]])
                        targets.append(node_to_idx[row[next_level]])
                        values.append(float(row["Weight"]))

                        if current_level == "TF":
                            all_connections.append(row[current_level])
                        else:
                            mask = flow_data[next_level] == row[next_level]
                            if current_level != "HAR":
                                mask &= flow_data[current_level] == row[current_level]
                            path_tfs = flow_data[mask]["TF"].unique()
                            if len(path_tfs) > 0:
                                tf_weights = {
                                    tf: flow_data[(flow_data["TF"] == tf) & mask][
                                        "Weight"
                                    ].sum()
                                    for tf in path_tfs
                                }
                                dominant_tf = max(
                                    tf_weights.items(), key=lambda x: x[1]
                                )[0]
                                all_connections.append(dominant_tf)

        n_tfs = len(top_tfs)
        colors = self._generate_colors(n_tfs)
        tf_color_map = {tf: colors[i % len(colors)] for i, tf in enumerate(top_tfs)}

        link_colors = [tf_color_map[tf] for tf in all_connections]

        node_colors = []
        for level in levels_to_show:
            if level == "TF":
                node_colors.extend(colors[: len(level_nodes[level])])
            elif level == "HAR":
                for har_group in level_nodes["HAR"]:
                    tf = har_to_tf_map[har_group]
                    node_colors.append(tf_color_map[tf])
            else:
                for node in level_nodes[level]:
                    mask = flow_data[level] == node
                    path_tfs = flow_data[mask]["TF"].unique()
                    if len(path_tfs) > 0:
                        tf_weights = {
                            tf: flow_data[(flow_data["TF"] == tf) & mask][
                                "Weight"
                            ].sum()
                            for tf in path_tfs
                        }
                        dominant_tf = max(tf_weights.items(), key=lambda x: x[1])[0]
                        node_colors.append(tf_color_map[dominant_tf])
                    else:
                        node_colors.append("rgba(200, 200, 200, 0.8)")

        fig = go.Figure(
            data=[
                go.Sankey(
                    node=dict(
                        pad=15,
                        thickness=20,
                        line=dict(color="black", width=0.5),
                        label=nodes,
                        color=node_colors,
                    ),
                    link=dict(
                        source=sources,
                        target=targets,
                        value=values,
                        color=link_colors,
                    ),
                )
            ]
        )

        title_parts = ["Enhanced Network Wiring"]
        if specific_celltypes:
            title_parts.append(f"Cell Types: {', '.join(specific_celltypes)}")
        if specific_regions:
            title_parts.append(f"Regions: {', '.join(specific_regions)}")
        if specific_stages:
            title_parts.append(f"Stages: {', '.join(specific_stages)}")

        title_text = " | ".join(title_parts)

        fig.update_layout(
            title=dict(
                text=title_text,
                x=0.1,
                y=0.95,
                xanchor="left",
                yanchor="top",
                font=dict(size=14, color="black", family="Arial"),
            ),
            font=dict(size=10, color="black", family="Arial"),
            height=height,
            width=width,
        )

        return fig, "Success"

    def get_node_subnetwork(self, node_name: str, node_type: str, max_depth: int = 1):
        """获取节点子网络，优化性能"""
        try:
            # 数据过滤优化
            if node_type == "tf":
                subnetwork_data = self.network_data[
                    (self.network_data["TF"] == node_name)
                    & (self.network_data["Weight"] > 0)
                ].copy()
            elif node_type == "target":
                subnetwork_data = self.network_data[
                    (self.network_data["Target"] == node_name)
                    & (self.network_data["Weight"] > 0)
                ].copy()
            else:
                return None

            if len(subnetwork_data) == 0:
                return None

            # 限制子网络大小以提高性能
            if len(subnetwork_data) > 100:
                subnetwork_data = subnetwork_data.nlargest(100, "Weight")

            G = nx.DiGraph()

            # 批量添加边，提高性能
            edges_data = []
            for _, row in subnetwork_data.iterrows():
                edges_data.append((
                    row["TF"],
                    row["Target"],
                    {
                        "weight": row["Weight"],
                        "celltype": row["CellType"],
                        "stage": row["Stage"],
                        "region": row["Region"],
                    }
                ))
            
            G.add_edges_from(edges_data)

            # 简化子图提取
            if node_name in G.nodes():
                try:
                    subgraph = nx.ego_graph(G, node_name, radius=max_depth, undirected=False)
                except:
                    subgraph = G
            else:
                subgraph = G

            # 限制节点数量
            if len(subgraph.nodes()) > 50:
                # 选择权重最高的连接
                node_weights = {}
                for node in subgraph.nodes():
                    total_weight = sum(data.get('weight', 0) for source, target, data in subgraph.edges(data=True) 
                                     if source == node or target == node)
                    node_weights[node] = total_weight
                
                top_nodes = sorted(node_weights.items(), key=lambda x: x[1], reverse=True)[:50]
                subgraph = subgraph.subgraph([node for node, _ in top_nodes])

            # 使用更快的布局算法
            try:
                pos = nx.spring_layout(subgraph, k=2, iterations=30, seed=42)
            except:
                pos = nx.random_layout(subgraph, seed=42)

            # 创建边轨迹
            edge_x = []
            edge_y = []
            edge_info = []

            for edge in subgraph.edges():
                x0, y0 = pos[edge[0]]
                x1, y1 = pos[edge[1]]
                edge_x.extend([x0, x1, None])
                edge_y.extend([y0, y1, None])

                edge_data = subgraph[edge[0]][edge[1]]
                edge_info.append({
                    "source": edge[0],
                    "target": edge[1],
                    "weight": edge_data.get("weight", 0),
                    "celltype": edge_data.get("celltype", "Unknown"),
                    "stage": edge_data.get("stage", "Unknown"),
                    "region": edge_data.get("region", "Unknown"),
                })

            # 创建节点轨迹
            node_x = []
            node_y = []
            node_text = []
            node_info = []

            for node in subgraph.nodes():
                x, y = pos[node]
                node_x.append(x)
                node_y.append(y)
                node_text.append(node)

                # 简化节点类型判断
                is_tf = any(edge[0] == node for edge in subgraph.edges())
                node_info.append({
                    "name": node,
                    "type": "TF" if is_tf else "Target",
                    "connections": subgraph.degree(node),
                })

            # 创建图表
            edge_trace = go.Scatter(
                x=edge_x,
                y=edge_y,
                line=dict(width=1.5, color="#888"),
                hoverinfo="none",
                mode="lines",
            )

            node_trace = go.Scatter(
                x=node_x,
                y=node_y,
                mode="markers+text",
                hoverinfo="text",
                text=node_text,
                textposition="middle center",
                marker=dict(
                    showscale=True,
                    colorscale="Viridis",
                    reversescale=False,
                    color=[],
                    size=15,
                    colorbar=dict(thickness=10, xanchor="left", titleside="right"),
                    line=dict(width=1, color="black"),
                ),
            )

            # 计算节点连接度
            node_adjacencies = []
            node_text_info = []
            for node in subgraph.nodes():
                degree = subgraph.degree(node)
                node_adjacencies.append(degree)
                node_text_info.append(f"{node}<br>Connections: {degree}")

            node_trace.marker.color = node_adjacencies
            node_trace.text = node_text_info

            fig = go.Figure(
                data=[edge_trace, node_trace],
                layout=go.Layout(
                    title=f"Subnetwork for {node_name} ({len(subgraph.nodes())} nodes, {len(subgraph.edges())} edges)",
                    titlefont_size=14,
                    showlegend=False,
                    hovermode="closest",
                    margin=dict(b=20, l=5, r=5, t=40),
                    annotations=[
                        dict(
                            text="Drag to explore • Zoom to see details",
                            showarrow=False,
                            xref="paper",
                            yref="paper",
                            x=0.005,
                            y=-0.002,
                            xanchor="left",
                            yanchor="bottom",
                            font=dict(color="gray", size=10),
                        )
                    ],
                    xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                    yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                ),
            )

            return fig, edge_info, node_info

        except Exception as e:
            print(f"Error in get_node_subnetwork: {str(e)}")
            return None

    def _generate_colors(self, n):
        if n <= 0:
            return []

        colors = []
        base_colors = [
            "rgb(25, 50, 200)",
            "rgb(205, 50, 65)",
            "rgb(255, 127, 14)",
            "rgb(44, 160, 44)",
            "rgb(25, 190, 200)",
            "rgb(148, 103, 189)",
            "rgb(227, 119, 194)",
            "rgb(188, 189, 34)",
        ]

        while len(colors) < n:
            alpha = 0.8 - (0.3 * (len(colors) // len(base_colors)))
            color_idx = len(colors) % len(base_colors)
            base_color = (
                base_colors[color_idx]
                .replace("rgb", "rgba")
                .replace(")", f", {alpha})")
            )
            colors.append(base_color)

        return colors

    def export_network_data(
        self,
        show_har: bool = True,
        top_n_tfs: int = 10,
        top_n_targets: int = 20,
        specific_tfs: List[str] = None,
        specific_targets: List[str] = None,
        specific_stages: List[str] = None,
        specific_regions: List[str] = None,
        specific_celltypes: List[str] = None,
        output_format: str = "csv",
    ):
        """
        Export filtered network data to CSV or JSON format

        Parameters:
        -----------
        output_format : str
            Export format ('csv' or 'json')

        Returns:
        --------
        str : Exported data as string
        """
        flow_data = self.network_data[
            (self.network_data["Weight"] > 0)
            & (self.network_data["TF"] != "Unknown_TF")
            & (self.network_data["Target"] != "Unknown_Target")
        ].copy()

        if specific_celltypes:
            flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]
        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]
        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]

        if specific_tfs:
            available_tfs = set(flow_data["TF"].unique())
            top_tfs = [tf for tf in specific_tfs if tf in available_tfs]
        else:
            top_tfs = self.get_top_nodes(
                "tf",
                top_n_tfs,
                {
                    "CellType": specific_celltypes,
                    "Region": specific_regions,
                    "Stage": specific_stages,
                },
            )

        if specific_targets:
            available_targets = set(flow_data["Target"].unique())
            top_targets = [
                target for target in specific_targets if target in available_targets
            ]
        else:
            flow_data_filtered = flow_data[flow_data["TF"].isin(top_tfs)]
            top_targets = self.get_top_nodes(
                "target",
                top_n_targets,
                {
                    "CellType": specific_celltypes,
                    "Region": specific_regions,
                    "Stage": specific_stages,
                },
            )

        export_data = flow_data[
            flow_data["TF"].isin(top_tfs) & flow_data["Target"].isin(top_targets)
        ].copy()

        if show_har and self.har_tf_data is not None:
            har_mapping = self.har_tf_data[self.har_tf_data["TF"].isin(top_tfs)]
            export_data = export_data.merge(har_mapping, on="TF", how="left")

        if output_format == "csv":
            output = io.StringIO()
            export_data.to_csv(output, index=False)
            return output.getvalue()
        elif output_format == "json":
            return export_data.to_json(orient="records", indent=2)
        else:
            raise ValueError("output_format must be 'csv' or 'json'")

    def get_network_statistics(
        self,
        specific_tfs: List[str] = None,
        specific_targets: List[str] = None,
        specific_stages: List[str] = None,
        specific_regions: List[str] = None,
        specific_celltypes: List[str] = None,
    ):
        """
        Get comprehensive network statistics

        Returns:
        --------
        dict : Network statistics
        """
        flow_data = self.network_data[
            (self.network_data["Weight"] > 0)
            & (self.network_data["TF"] != "Unknown_TF")
            & (self.network_data["Target"] != "Unknown_Target")
        ].copy()

        if specific_celltypes:
            flow_data = flow_data[flow_data["CellType"].isin(specific_celltypes)]
        if specific_regions:
            flow_data = flow_data[flow_data["Region"].isin(specific_regions)]
        if specific_stages:
            flow_data = flow_data[flow_data["Stage"].isin(specific_stages)]
        if specific_tfs:
            flow_data = flow_data[flow_data["TF"].isin(specific_tfs)]
        if specific_targets:
            flow_data = flow_data[flow_data["Target"].isin(specific_targets)]

        stats = {
            "total_connections": len(flow_data),
            "unique_tfs": flow_data["TF"].nunique(),
            "unique_targets": flow_data["Target"].nunique(),
            "unique_celltypes": flow_data["CellType"].nunique(),
            "unique_regions": flow_data["Region"].nunique(),
            "unique_stages": flow_data["Stage"].nunique(),
            "total_weight": flow_data["Weight"].sum(),
            "avg_weight": flow_data["Weight"].mean(),
            "max_weight": flow_data["Weight"].max(),
            "min_weight": flow_data["Weight"].min(),
            "top_tfs_by_weight": flow_data.groupby("TF")["Weight"]
            .sum()
            .nlargest(10)
            .to_dict(),
            "top_targets_by_weight": flow_data.groupby("Target")["Weight"]
            .sum()
            .nlargest(10)
            .to_dict(),
            "connections_by_celltype": flow_data.groupby("CellType").size().to_dict(),
            "connections_by_region": flow_data.groupby("Region").size().to_dict(),
            "connections_by_stage": flow_data.groupby("Stage").size().to_dict(),
        }

        return stats


def create_dash_app(visualizer: EnhancedNetworkVisualizer):
    app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP])

    options = visualizer.get_available_options()

    app.layout = dbc.Container(
        [
            dbc.Row(
                [
                    dbc.Col(
                        [
                            html.H1(
                                "Enhanced Network Visualization",
                                className="text-center mb-4",
                            ),
                            html.Hr(),
                        ]
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Card(
                                [
                                    dbc.CardHeader("Filter Controls"),
                                    dbc.CardBody(
                                        [
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                "Show HAR Level:"
                                                            ),
                                                            dbc.Switch(
                                                                id="show-har-switch",
                                                                value=True,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Top N TFs:"),
                                                            dbc.Input(
                                                                id="top-n-tfs",
                                                                type="number",
                                                                value=10,
                                                                min=1,
                                                                max=50,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ]
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Top N Targets:"),
                                                            dbc.Input(
                                                                id="top-n-targets",
                                                                type="number",
                                                                value=20,
                                                                min=1,
                                                                max=100,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Chart Height:"),
                                                            dbc.Input(
                                                                id="chart-height",
                                                                type="number",
                                                                value=600,
                                                                min=300,
                                                                max=1000,
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ]
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Cell Types:"),
                                                            dcc.Dropdown(
                                                                id="celltype-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": ct,
                                                                        "value": ct,
                                                                    }
                                                                    for ct in options[
                                                                        "celltypes"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select cell types...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Regions:"),
                                                            dcc.Dropdown(
                                                                id="region-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": r,
                                                                        "value": r,
                                                                    }
                                                                    for r in options[
                                                                        "regions"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select regions...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Stages:"),
                                                            dcc.Dropdown(
                                                                id="stage-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": s,
                                                                        "value": s,
                                                                    }
                                                                    for s in options[
                                                                        "stages"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select stages...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Specific TFs:"),
                                                            dcc.Dropdown(
                                                                id="tf-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": tf,
                                                                        "value": tf,
                                                                    }
                                                                    for tf in options[
                                                                        "tfs"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select specific TFs...",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label(
                                                                "Specific Targets:"
                                                            ),
                                                            dcc.Dropdown(
                                                                id="target-dropdown",
                                                                options=[
                                                                    {
                                                                        "label": t,
                                                                        "value": t,
                                                                    }
                                                                    for t in options[
                                                                        "targets"
                                                                    ]
                                                                ],
                                                                multi=True,
                                                                placeholder="Select specific targets...",
                                                            ),
                                                        ],
                                                        width=12,
                                                    )
                                                ],
                                                className="mt-3",
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Button(
                                                                "Generate Sankey Diagram",
                                                                id="generate-button",
                                                                color="primary",
                                                                className="mt-3 w-100",
                                                            )
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Button(
                                                                "Export Data (CSV)",
                                                                id="export-csv-button",
                                                                color="success",
                                                                className="mt-3 w-100",
                                                            )
                                                        ],
                                                        width=6,
                                                    ),
                                                ]
                                            ),
                                        ]
                                    ),
                                ]
                            )
                        ],
                        width=4,
                    ),
                    dbc.Col(
                        [
                            dbc.Card(
                                [
                                    dbc.CardHeader("Network Visualization"),
                                    dbc.CardBody(
                                        [
                                            dcc.Graph(
                                                id="sankey-plot",
                                                style={"height": "600px"},
                                            ),
                                            html.Div(id="click-info", className="mt-3"),
                                            html.Div(
                                                id="download-container",
                                                className="mt-3",
                                            ),
                                        ]
                                    ),
                                ]
                            )
                        ],
                        width=8,
                    ),
                ],
                className="mb-4",
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Card(
                                [
                                    dbc.CardHeader("Node Subnetwork Explorer"),
                                    dbc.CardBody(
                                        [
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Node Name:"),
                                                            dbc.Input(
                                                                id="node-name-input",
                                                                placeholder="Enter node name...",
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                    dbc.Col(
                                                        [
                                                            dbc.Label("Node Type:"),
                                                            dbc.Select(
                                                                id="node-type-select",
                                                                options=[
                                                                    {
                                                                        "label": "TF",
                                                                        "value": "tf",
                                                                    },
                                                                    {
                                                                        "label": "Target",
                                                                        "value": "target",
                                                                    },
                                                                ],
                                                                value="tf",
                                                                className="mb-3",
                                                            ),
                                                        ],
                                                        width=6,
                                                    ),
                                                ]
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dbc.Button(
                                                                "Show Subnetwork",
                                                                id="subnetwork-button",
                                                                color="secondary",
                                                                className="w-100",
                                                            )
                                                        ]
                                                    )
                                                ]
                                            ),
                                            dbc.Row(
                                                [
                                                    dbc.Col(
                                                        [
                                                            dcc.Graph(
                                                                id="subnetwork-plot",
                                                                style={
                                                                    "height": "400px"
                                                                },
                                                            )
                                                        ]
                                                    )
                                                ],
                                                className="mt-3",
                                            ),
                                        ]
                                    ),
                                ]
                            )
                        ],
                        width=12,
                    )
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(
                        [
                            dbc.Card(
                                [
                                    dbc.CardHeader("Network Statistics"),
                                    dbc.CardBody([html.Div(id="network-stats")]),
                                ]
                            )
                        ],
                        width=12,
                    )
                ],
                className="mt-4",
            ),
        ],
        fluid=True,
    )

    @app.callback(
        [Output("sankey-plot", "figure"), Output("click-info", "children")],
        [Input("generate-button", "n_clicks")],
        [
            State("show-har-switch", "value"),
            State("top-n-tfs", "value"),
            State("top-n-targets", "value"),
            State("celltype-dropdown", "value"),
            State("region-dropdown", "value"),
            State("stage-dropdown", "value"),
            State("tf-dropdown", "value"),
            State("target-dropdown", "value"),
            State("chart-height", "value"),
        ],
    )
    def update_sankey_plot(
        n_clicks,
        show_har,
        top_n_tfs,
        top_n_targets,
        celltypes,
        regions,
        stages,
        tfs,
        targets,
        height,
    ):
        if n_clicks is None:
            return go.Figure(), ""

        fig, message = visualizer.create_sankey_diagram(
            show_har=show_har,
            top_n_tfs=top_n_tfs,
            top_n_targets=top_n_targets,
            specific_tfs=tfs,
            specific_targets=targets,
            specific_stages=stages,
            specific_regions=regions,
            specific_celltypes=celltypes,
            height=height,
        )

        if fig is None:
            return go.Figure(), dbc.Alert(message, color="warning")

        return fig, dbc.Alert("Sankey diagram generated successfully!", color="success")

    @app.callback(
        Output("download-container", "children"),
        [Input("export-csv-button", "n_clicks")],
        [
            State("show-har-switch", "value"),
            State("top-n-tfs", "value"),
            State("top-n-targets", "value"),
            State("celltype-dropdown", "value"),
            State("region-dropdown", "value"),
            State("stage-dropdown", "value"),
            State("tf-dropdown", "value"),
            State("target-dropdown", "value"),
        ],
    )
    def export_data(
        n_clicks,
        show_har,
        top_n_tfs,
        top_n_targets,
        celltypes,
        regions,
        stages,
        tfs,
        targets,
    ):
        if n_clicks is None:
            return ""

        try:
            csv_data = visualizer.export_network_data(
                show_har=show_har,
                top_n_tfs=top_n_tfs,
                top_n_targets=top_n_targets,
                specific_tfs=tfs,
                specific_targets=targets,
                specific_stages=stages,
                specific_regions=regions,
                specific_celltypes=celltypes,
                output_format="csv",
            )

            encoded_data = base64.b64encode(csv_data.encode()).decode()

            return html.Div(
                [
                    html.A(
                        "Download CSV File",
                        id="download-link",
                        download="network_data.csv",
                        href=f"data:text/csv;base64,{encoded_data}",
                        className="btn btn-success",
                    )
                ]
            )
        except Exception as e:
            return dbc.Alert(f"Export failed: {str(e)}", color="danger")

    @app.callback(
        [Output("subnetwork-plot", "figure"), Output("network-stats", "children")],
        [Input("subnetwork-button", "n_clicks")],
        [State("node-name-input", "value"), State("node-type-select", "value")],
    )
    def update_subnetwork(n_clicks, node_name, node_type):
        if n_clicks is None or not node_name:
            return go.Figure(), ""

        try:
            result = visualizer.get_node_subnetwork(node_name, node_type)
            if result is None:
                return go.Figure(), dbc.Alert(
                    f"No subnetwork found for node '{node_name}' of type '{node_type}'.", 
                    color="warning"
                )

            fig, edge_info, node_info = result

            if not node_info:
                return fig, dbc.Alert("No node information available.", color="warning")

            # 创建统计表格
            stats_table = dash_table.DataTable(
                data=node_info,
                columns=[{"name": i, "id": i} for i in node_info[0].keys()],
                style_cell={"textAlign": "left", "fontSize": 12},
                style_header={
                    "backgroundColor": "rgb(230, 230, 230)",
                    "fontWeight": "bold",
                },
                style_data_conditional=[
                    {
                        'if': {'row_index': 'odd'},
                        'backgroundColor': 'rgb(248, 248, 248)'
                    }
                ],
                page_size=10,
            )

            return fig, stats_table

        except Exception as e:
            return go.Figure(), dbc.Alert(
                f"Error generating subnetwork: {str(e)}", 
                color="danger"
            )

    return app


if __name__ == "__main__":
    try:
        data_path = "data/networks/csv/network_data.csv"
        har_tf_path = "data/networks/csv/tf_har_combinations.csv"

        print(f"{get_colored_time()} Initializing Enhanced Network Visualizer...")
        visualizer = EnhancedNetworkVisualizer(data_path, har_tf_path)

        print(f"{get_colored_time()} Creating Dash application...")
        app = create_dash_app(visualizer)

        print(f"{get_colored_time()} Starting web server...")
        print(
            f"{get_colored_time()} Open your browser and go to: http://127.0.0.1:8000"
        )

        app.run(debug=True, host="127.0.0.1", port=8000)

    except Exception as e:
        print(f"{get_colored_time()} {Fore.RED}Error: {str(e)}{Style.RESET_ALL}")
        import traceback

        traceback.print_exc()
        sys.exit(1)
