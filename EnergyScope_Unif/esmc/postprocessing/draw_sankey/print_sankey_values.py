import sys, getopt
from pathlib import Path
import pandas as pd
import plotly.graph_objects as go
import numpy as np

def read_and_process_data(filepath, num_decimals):
    filepath = Path(filepath)
    df = pd.read_csv(str(filepath / "input2sankey_Total.csv"))

    source_sum = df.groupby('source')['realValue'].sum().to_dict()
    target_sum = df.groupby('target')['realValue'].sum().to_dict()

    df['source'] = df['source'].apply(lambda x: f"{x} ({round(source_sum[x], num_decimals)})")
    df['target'] = df['target'].apply(lambda x: f"{x} ({round(source_sum[x], num_decimals)})" if x in source_sum else x)
    df['target'] = df.apply(lambda row: f"{row['target']} ({round(target_sum[row['target']], num_decimals)})" if row['target'] not in df['source'].unique() else row['target'], axis=1)

    labels = pd.concat([df['source'], df['target']]).unique()
    np.random.shuffle(labels)  # Shuffle labels
    label_dict = {label: idx for idx, label in enumerate(labels)}
    df['source'] = df['source'].map(label_dict)
    df['target'] = df['target'].map(label_dict)

    return df, labels

def generate_sankey(df, labels, link_transparency, node_color, filepath, auto_open=True):
    filepath = Path(filepath)
    color_dict = {i: f'rgba({color.strip("#")},0.5)' for i, color in df['layerColor'].items()}
    unit_dict = {label: unit for label, unit in zip(df['source'], df['layerUnit'])}
    units = [unit_dict.get(label, "") for label in labels]

    fig = go.Figure(data=[go.Sankey(
        node = dict(
          pad = 15,
          thickness = 20,
          line = dict(color = "black", width = 0.5),
          label = labels,
          color = node_color
        ),
        link = dict(
          source = df['source'],
          target = df['target'],
          value = df['realValue'],
          label = units,
          color = df['layerColor'].apply(lambda x: f'rgba({",".join(map(str, list(int(x[i:i+2], 16) for i in (1, 3, 5))))},{link_transparency})')
        )
    )])

    fig.update_layout(title_text="Sankey Diagram: Bolivia", font_size=18)
    fig.write_html(str(filepath / "generated_sankey_Total.html"), auto_open=auto_open)
    #fig.show()

if __name__ == "__main__":
    main()
