import matplotlib.pyplot as plt


def plot_donqin_airfoil(df_all, output_path, show=False):
    # Plot
    fig = plt.figure(figsize=(16, 6.5))
    ax = fig.add_subplot(111)

    pieces = df_all['piece_id'].unique()
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    for piece_i, color_i in zip(pieces, colors):
        df_i = df_all.copy()[df_all['piece_id'] == piece_i]

        db_ini = df_i.db.iloc[0]
        db_end = df_i.db.iloc[-1]
        piece_label = f'{piece_i} - db{db_ini} to db{db_end}'
        ini_text = f'ini_{piece_i} (db{db_ini})'
        end_text = f'end_{piece_i} (db{db_end})'
        ax.plot(df_i.x, df_i.y, '-o', color=color_i, label=piece_label)
        ax.text(df_i.x.iloc[0], df_i.y.iloc[0]+0.025, ini_text, color=color_i)
        ax.text(df_i.x.iloc[-1], df_i.y.iloc[-1]+0.0, end_text, color=color_i)

    ax.grid()
    ax.legend()
    ax.set_xlabel('X')
    ax.set_xlabel('Y')
    ax.set_ylim(-0.15, 0.22)  # same as LEI
    ax.set_xlim(-0.05, 1.05)  # same as LEI
    ax.set_aspect('equal', adjustable='box')
    fig.suptitle('Donqin duct')

    # Write
    fig.savefig(output_path, dpi=600)
    if show:
        plt.show()

