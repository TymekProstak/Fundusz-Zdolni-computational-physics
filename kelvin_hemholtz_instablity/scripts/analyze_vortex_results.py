#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
import shutil
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image


@dataclass
class StateFile:
    delta_text: str
    delta_value: float
    step: int
    time: float
    path: Path


# =============================================================================
#                               WCZYTYWANIE DANYCH
# =============================================================================

def parse_result_filename(path: Path, prefix: str):
    pattern = rf"{prefix}_(.+)_(\d+)\.txt$"
    match = re.match(pattern, path.name)

    if match is None:
        return None

    delta_text = match.group(1)
    step = int(match.group(2))

    try:
        delta_value = float(delta_text)
    except ValueError:
        return None

    return delta_text, delta_value, step


def calculate_time_from_step(step: int, delta_value: float, L: float, u1: float, u2: float):
    circulation = (u1 - u2) * L

    if abs(circulation) < 1e-14:
        return float(step)

    l_docelowe = 0.05 * delta_value
    dt = 0.250 * L * l_docelowe / abs(circulation)

    return step * dt


def find_state_files(input_dir: Path, L: float, u1: float, u2: float):
    files = []

    for path in input_dir.glob("wyniki_*.txt"):
        parsed = parse_result_filename(path, "wyniki")

        if parsed is None:
            continue

        delta_text, delta_value, step = parsed
        time = calculate_time_from_step(step, delta_value, L, u1, u2)

        files.append(
            StateFile(
                delta_text=delta_text,
                delta_value=delta_value,
                step=step,
                time=time,
                path=path,
            )
        )

    files.sort(key=lambda item: (item.delta_value, item.step))

    return files


def group_by_delta(files):
    groups = {}

    for item in files:
        if item.delta_text not in groups:
            groups[item.delta_text] = []

        groups[item.delta_text].append(item)

    return groups


def load_half_file(path: Path):
    values = np.loadtxt(path, dtype=float)

    if values.ndim != 1:
        values = values.reshape(-1)

    if len(values) % 2 != 0:
        raise ValueError(f"Nieparzysta liczba wartosci w pliku: {path}")

    n = len(values) // 2

    first = values[:n]
    second = values[n:]

    return first, second


def load_sheet(path: Path):
    x, y = load_half_file(path)
    return x, y


def load_velocity_file(input_dir: Path, delta_text: str, step: int):
    path = input_dir / f"predkosc_{delta_text}_{step}.txt"

    if not path.exists():
        return None

    ux, uy = load_half_file(path)
    return ux, uy


# =============================================================================
#                               SCIEZKI
# =============================================================================

def resolve_paths(input_arg: str, output_arg: str):
    input_dir = Path(input_arg)
    output_dir = Path(output_arg)

    # jezeli skrypt jest odpalony z folderu scripts/
    # i nie znajduje results/, to sprawdza ../results/
    if not input_dir.exists():
        parent_input = Path("..") / input_arg

        if parent_input.exists():
            input_dir = parent_input

            if output_arg == "figures":
                output_dir = Path("..") / output_arg

    return input_dir, output_dir


# =============================================================================
#                               WYBOR STANOW
# =============================================================================

def choose_deltas(groups, requested_deltas):
    if requested_deltas:
        selected = []

        for delta in requested_deltas:
            if delta in groups:
                selected.append(delta)
            else:
                print(f"[warning] Brak danych dla delta={delta}")

        return selected

    return sorted(groups.keys(), key=lambda d: float(d))


def choose_states_for_plots(states, requested_times, snapshot_count):
    if not states:
        return []

    selected = []

    if requested_times:
        for target_time in requested_times:
            closest = min(states, key=lambda item: abs(item.time - target_time))
            selected.append(closest)
    else:
        if len(states) <= snapshot_count:
            selected = states
        else:
            indices = np.linspace(0, len(states) - 1, snapshot_count)
            indices = np.round(indices).astype(int)
            selected = [states[i] for i in indices]

    unique = {}

    for item in selected:
        unique[item.step] = item

    selected = list(unique.values())
    selected.sort(key=lambda item: item.step)

    return selected


def downsample_states(states, max_count):
    if len(states) <= max_count:
        return states

    indices = np.linspace(0, len(states) - 1, max_count)
    indices = np.round(indices).astype(int)

    return [states[i] for i in indices]


def auto_compare_time(groups, selected_deltas):
    # wspolny maksymalny czas, zeby porownanie mialo sens dla wszystkich delt
    last_times = []

    for delta_text in selected_deltas:
        states = groups[delta_text]

        if states:
            last_times.append(states[-1].time)

    if not last_times:
        return 0.0

    return min(last_times)


def print_selected_time_info(groups, selected_deltas):
    for delta_text in selected_deltas:
        states = groups[delta_text]

        if not states:
            continue

        print(
            f"[delta={delta_text}] zakres czasu: "
            f"{states[0].time:.6g} -> {states[-1].time:.6g}, "
            f"kroki: {states[0].step} -> {states[-1].step}"
        )


# =============================================================================
#                               ANALIZA FOURIEROWSKA
# =============================================================================

def spectrum_from_velocity(ux, uy):
    ux = np.asarray(ux, dtype=float)
    uy = np.asarray(uy, dtype=float)

    ux = ux - np.mean(ux)
    uy = uy - np.mean(uy)

    n = len(ux)

    fft_x = np.fft.rfft(ux)
    fft_y = np.fft.rfft(uy)

    energy = (np.abs(fft_x) ** 2 + np.abs(fft_y) ** 2) / (n * n)
    modes = np.arange(len(energy))

    return modes, energy


def spectrum_from_geometry(y):
    y = np.asarray(y, dtype=float)
    y = y - np.mean(y)

    n = len(y)

    fft_y = np.fft.rfft(y)
    energy = (np.abs(fft_y) ** 2) / (n * n)
    modes = np.arange(len(energy))

    return modes, energy


def calculate_flow(uy, L):
    uy = np.asarray(uy, dtype=float)
    x = np.linspace(0.0, L, len(uy))

    flow_abs = np.trapz(np.abs(uy), x)
    flow_signed = np.trapz(uy, x)

    return flow_abs, flow_signed


# =============================================================================
#                               WYKRESY
# =============================================================================

def get_sheet_limits(states):
    all_x = []
    all_y = []

    for state in states:
        x, y = load_sheet(state.path)
        all_x.append(x)
        all_y.append(y)

    all_x = np.concatenate(all_x)
    all_y = np.concatenate(all_y)

    x_span = np.max(all_x) - np.min(all_x)
    y_span = np.max(all_y) - np.min(all_y)

    if x_span < 1e-14:
        x_span = 1.0

    if y_span < 1e-14:
        y_span = 1.0

    x_min = np.min(all_x) - 0.05 * x_span
    x_max = np.max(all_x) + 0.05 * x_span

    y_min = np.min(all_y) - 0.20 * y_span
    y_max = np.max(all_y) + 0.20 * y_span

    return x_min, x_max, y_min, y_max


def plot_sheet_snapshots(delta_text, selected_states, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)

    n = len(selected_states)

    fig, axes = plt.subplots(
        n,
        1,
        figsize=(10, 2.8 * n),
        sharex=True
    )

    if n == 1:
        axes = [axes]

    x_min, x_max, y_min, y_max = get_sheet_limits(selected_states)

    for ax, state in zip(axes, selected_states):
        x, y = load_sheet(state.path)

        marker_step = max(1, len(x) // 400)

        ax.plot(x, y, linewidth=1.2)
        ax.scatter(x[::marker_step], y[::marker_step], s=8)

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        ax.grid(True, alpha=0.3)
        ax.set_ylabel("y")

        ax.set_title(
            f"delta={delta_text}, step={state.step}, t={state.time:.4g}",
            fontsize=10,
            pad=8
        )

    axes[-1].set_xlabel("x")

    fig.suptitle(
        f"Vortex sheet evolution, delta={delta_text}",
        fontsize=15,
        y=0.995
    )

    fig.tight_layout(rect=[0.0, 0.0, 1.0, 0.955])

    path = output_dir / f"sheet_snapshots_delta_{delta_text}.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)

    return path


def plot_spectrum_snapshots(delta_text, selected_states, input_dir, output_dir):
    output_dir.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 5))

    used_any = False

    for state in selected_states:
        velocity = load_velocity_file(input_dir, delta_text, state.step)

        if velocity is not None:
            ux, uy = velocity
            modes, energy = spectrum_from_velocity(ux, uy)
            label = f"t={state.time:.4g}, step={state.step}"
            used_any = True
        else:
            _, y = load_sheet(state.path)
            modes, energy = spectrum_from_geometry(y)
            label = f"t={state.time:.4g}, step={state.step} (geometry)"

        if len(modes) > 1:
            ax.semilogy(modes[1:], energy[1:] + 1e-30, linewidth=1.5, label=label)

    ax.set_xlabel("Fourier mode number")
    ax.set_ylabel("spectral energy")
    ax.set_title(f"Fourier spectrum, delta={delta_text}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=8)

    if not used_any:
        ax.text(
            0.02,
            0.02,
            "No predkosc_*.txt files found - spectrum of geometry y(x) is shown.",
            transform=ax.transAxes,
            fontsize=9,
        )

    fig.tight_layout()

    path = output_dir / f"spectrum_snapshots_delta_{delta_text}.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)

    return path


def plot_flow(delta_text, states, input_dir, output_dir, L):
    output_dir.mkdir(parents=True, exist_ok=True)

    times = []
    flow_abs_values = []
    flow_signed_values = []

    for state in states:
        velocity = load_velocity_file(input_dir, delta_text, state.step)

        if velocity is None:
            continue

        _, uy = velocity
        flow_abs, flow_signed = calculate_flow(uy, L)

        times.append(state.time)
        flow_abs_values.append(flow_abs)
        flow_signed_values.append(flow_signed)

    if not times:
        return None

    fig, ax = plt.subplots(figsize=(9, 5))

    ax.plot(times, flow_abs_values, linewidth=1.8, label="integral |Uy| along y=0")
    ax.plot(times, flow_signed_values, linewidth=1.2, label="integral Uy along y=0")

    ax.set_xlabel("time")
    ax.set_ylabel("flow")
    ax.set_title(f"Flow through y=0, delta={delta_text}")
    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()

    path = output_dir / f"flow_delta_{delta_text}.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)

    return path


def plot_delta_comparison(groups, selected_deltas, input_dir, output_dir, L, target_time):
    output_dir.mkdir(parents=True, exist_ok=True)

    fig_sheet, ax_sheet = plt.subplots(figsize=(10, 5))
    fig_spectrum, ax_spectrum = plt.subplots(figsize=(9, 5))
    fig_flow, ax_flow = plt.subplots(figsize=(9, 5))

    spectrum_any = False
    flow_any = False

    actual_times = []

    for delta_text in selected_deltas:
        states = groups[delta_text]

        state = min(states, key=lambda item: abs(item.time - target_time))
        actual_times.append(state.time)

        x, y = load_sheet(state.path)

        ax_sheet.plot(
            x,
            y,
            linewidth=1.3,
            label=f"delta={delta_text}, t={state.time:.4g}",
        )

        velocity = load_velocity_file(input_dir, delta_text, state.step)

        if velocity is not None:
            ux, uy = velocity
            modes, energy = spectrum_from_velocity(ux, uy)

            if len(modes) > 1:
                ax_spectrum.semilogy(
                    modes[1:],
                    energy[1:] + 1e-30,
                    linewidth=1.4,
                    label=f"delta={delta_text}, t={state.time:.4g}",
                )
                spectrum_any = True

        times = []
        flow_values = []

        for item in states:
            velocity_item = load_velocity_file(input_dir, delta_text, item.step)

            if velocity_item is None:
                continue

            _, uy_item = velocity_item
            flow_abs, _ = calculate_flow(uy_item, L)

            times.append(item.time)
            flow_values.append(flow_abs)

        if times:
            ax_flow.plot(times, flow_values, linewidth=1.5, label=f"delta={delta_text}")
            flow_any = True

    actual_time = min(actual_times) if actual_times else target_time

    ax_sheet.set_xlabel("x")
    ax_sheet.set_ylabel("y")
    ax_sheet.set_title(f"Vortex sheet comparison for t≈{actual_time:.4g}")
    ax_sheet.grid(True, alpha=0.3)
    ax_sheet.legend(fontsize=8)
    fig_sheet.tight_layout()

    sheet_path = output_dir / f"comparison_sheet_t_{actual_time:g}.png"
    fig_sheet.savefig(sheet_path, dpi=160)
    plt.close(fig_sheet)

    spectrum_path = None

    if spectrum_any:
        ax_spectrum.set_xlabel("Fourier mode number")
        ax_spectrum.set_ylabel("spectral energy")
        ax_spectrum.set_title(f"Fourier spectrum comparison for t≈{actual_time:.4g}")
        ax_spectrum.grid(True, alpha=0.3)
        ax_spectrum.legend(fontsize=8)
        fig_spectrum.tight_layout()

        spectrum_path = output_dir / f"comparison_spectrum_t_{actual_time:g}.png"
        fig_spectrum.savefig(spectrum_path, dpi=160)

    plt.close(fig_spectrum)

    flow_path = None

    if flow_any:
        ax_flow.set_xlabel("time")
        ax_flow.set_ylabel("integral |Uy| along y=0")
        ax_flow.set_title("Flow through y=0 comparison")
        ax_flow.grid(True, alpha=0.3)
        ax_flow.legend(fontsize=8)
        fig_flow.tight_layout()

        flow_path = output_dir / "comparison_flow_delta.png"
        fig_flow.savefig(flow_path, dpi=160)

    plt.close(fig_flow)

    return sheet_path, spectrum_path, flow_path


# =============================================================================
#                               ANIMACJA
# =============================================================================

def make_animation_for_delta(delta_text, states, input_dir, output_dir, fps, max_frames):
    output_dir.mkdir(parents=True, exist_ok=True)

    selected_states = downsample_states(states, max_frames)

    frame_dir = output_dir / f"_frames_delta_{delta_text}"

    if frame_dir.exists():
        shutil.rmtree(frame_dir)

    frame_dir.mkdir(parents=True, exist_ok=True)

    x_min, x_max, y_min, y_max = get_sheet_limits(selected_states)

    frame_paths = []

    for frame_id, state in enumerate(selected_states):
        x, y = load_sheet(state.path)

        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))

        ax_sheet = axes[0]
        ax_spec = axes[1]

        marker_step = max(1, len(x) // 500)

        ax_sheet.plot(x, y, linewidth=1.2)
        ax_sheet.scatter(x[::marker_step], y[::marker_step], s=6)

        ax_sheet.set_xlim(x_min, x_max)
        ax_sheet.set_ylim(y_min, y_max)

        ax_sheet.set_xlabel("x")
        ax_sheet.set_ylabel("y")
        ax_sheet.grid(True, alpha=0.3)
        ax_sheet.set_title(f"Vortex sheet\nstep={state.step}, t={state.time:.4g}")

        velocity = load_velocity_file(input_dir, delta_text, state.step)

        if velocity is not None:
            ux, uy = velocity
            modes, energy = spectrum_from_velocity(ux, uy)
            title = "Velocity spectrum"
        else:
            modes, energy = spectrum_from_geometry(y)
            title = "Geometry spectrum y(x)"

        if len(modes) > 1:
            ax_spec.semilogy(modes[1:], energy[1:] + 1e-30, linewidth=1.4)

        ax_spec.set_xlabel("mode number")
        ax_spec.set_ylabel("energy")
        ax_spec.grid(True, alpha=0.3)
        ax_spec.set_title(title)

        fig.suptitle(f"delta={delta_text}", fontsize=14)
        fig.tight_layout()

        frame_path = frame_dir / f"frame_{frame_id:05d}.png"
        fig.savefig(frame_path, dpi=130)
        plt.close(fig)

        frame_paths.append(frame_path)

        print(f"[animation delta={delta_text}] frame {frame_id + 1}/{len(selected_states)}")

    gif_path = output_dir / f"vortex_delta_{delta_text}.gif"

    images = [Image.open(path) for path in frame_paths]
    duration = int(1000 / fps)

    images[0].save(
        gif_path,
        save_all=True,
        append_images=images[1:],
        duration=duration,
        loop=0,
    )

    for image in images:
        image.close()

    # klatki sa tymczasowe, gif zostaje
    shutil.rmtree(frame_dir)

    return gif_path


# =============================================================================
#                               MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Analiza wynikow symulacji powierzchni wirowej."
    )

    parser.add_argument(
        "--input",
        type=str,
        default="results",
        help="Folder z plikami wyniki_*.txt oraz predkosc_*.txt.",
    )

    parser.add_argument(
        "--output",
        type=str,
        default="figures",
        help="Folder wyjsciowy na wykresy i animacje.",
    )

    parser.add_argument(
        "--L",
        type=float,
        default=1.0,
        help="Dlugosc poczatkowej powierzchni wirowej.",
    )

    parser.add_argument(
        "--u1",
        type=float,
        default=2.0,
        help="Predkosc po jednej stronie powierzchni.",
    )

    parser.add_argument(
        "--u2",
        type=float,
        default=1.0,
        help="Predkosc po drugiej stronie powierzchni.",
    )

    parser.add_argument(
        "--deltas",
        nargs="*",
        default=None,
        help="Lista wybranych delt, np. --deltas 0.5 0.25 0.05.",
    )

    parser.add_argument(
        "--times",
        nargs="*",
        type=float,
        default=None,
        help="Wybrane chwile czasu do wykresow, np. --times 0 1 2 4.5. Jezeli nie podano, chwile sa dobierane automatycznie.",
    )

    parser.add_argument(
        "--snapshot-count",
        type=int,
        default=5,
        help="Liczba automatycznie wybranych chwil czasu, jezeli nie podano --times.",
    )

    parser.add_argument(
        "--compare-time",
        type=float,
        default=None,
        help="Czas do porownania roznych delt. Jezeli nie podano, wybierany jest ostatni wspolny dostepny czas.",
    )

    parser.add_argument(
        "--fps",
        type=int,
        default=12,
        help="Liczba klatek na sekunde w GIF.",
    )

    parser.add_argument(
        "--max-animation-frames",
        type=int,
        default=140,
        help="Maksymalna liczba klatek animacji dla jednej delty.",
    )

    args = parser.parse_args()

    input_dir, output_dir = resolve_paths(args.input, args.output)

    plots_dir = output_dir / "plots"
    animations_dir = output_dir / "animations"
    comparisons_dir = output_dir / "comparisons"

    print(f"Input folder:  {input_dir}")
    print(f"Output folder: {output_dir}")

    files = find_state_files(input_dir, args.L, args.u1, args.u2)

    if not files:
        raise RuntimeError(
            f"Nie znaleziono plikow wyniki_*.txt w folderze: {input_dir}"
        )

    groups = group_by_delta(files)
    selected_deltas = choose_deltas(groups, args.deltas)

    if not selected_deltas:
        raise RuntimeError("Nie wybrano zadnej delty do analizy.")

    print("Znalezione delty:", ", ".join(selected_deltas))
    print_selected_time_info(groups, selected_deltas)

    if args.compare_time is None:
        args.compare_time = auto_compare_time(groups, selected_deltas)
        print(f"Czas porownania dobrany automatycznie: t = {args.compare_time:.6g}")
    else:
        print(f"Czas porownania podany recznie: t = {args.compare_time:.6g}")

    generated_files = []

    for delta_text in selected_deltas:
        states = groups[delta_text]

        selected_states = choose_states_for_plots(
            states=states,
            requested_times=args.times,
            snapshot_count=args.snapshot_count,
        )

        print(f"\n[delta={delta_text}] liczba stanow: {len(states)}")
        print("[delta={}] wybrane kroki: {}".format(
            delta_text,
            ", ".join(str(state.step) for state in selected_states),
        ))

        print("[delta={}] wybrane czasy: {}".format(
            delta_text,
            ", ".join(f"{state.time:.6g}" for state in selected_states),
        ))

        path = plot_sheet_snapshots(delta_text, selected_states, plots_dir)
        generated_files.append(path)

        path = plot_spectrum_snapshots(delta_text, selected_states, input_dir, plots_dir)
        generated_files.append(path)

        path = plot_flow(delta_text, states, input_dir, plots_dir, args.L)

        if path is not None:
            generated_files.append(path)
        else:
            print(f"[warning] Brak predkosc_*.txt dla delta={delta_text}; pomijam przeplyw.")

        path = make_animation_for_delta(
            delta_text=delta_text,
            states=states,
            input_dir=input_dir,
            output_dir=animations_dir,
            fps=args.fps,
            max_frames=args.max_animation_frames,
        )
        generated_files.append(path)

    comparison_paths = plot_delta_comparison(
        groups=groups,
        selected_deltas=selected_deltas,
        input_dir=input_dir,
        output_dir=comparisons_dir,
        L=args.L,
        target_time=args.compare_time,
    )

    for path in comparison_paths:
        if path is not None:
            generated_files.append(path)

    print("\nWygenerowane pliki:")

    for path in generated_files:
        print(f"  {path}")


if __name__ == "__main__":
    main()
