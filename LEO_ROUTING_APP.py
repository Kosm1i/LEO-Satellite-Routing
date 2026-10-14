from skyfield.api import load, wgs84
import numpy as np
import heapq
import time
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import tkinter as tk
from tkinter import scrolledtext, messagebox, ttk
import threading

MIASTA = {
    "Poznań": (52.4064, 16.9252),
    "Warszawa": (52.2297, 21.0122),
    "Nowy Jork": (40.7128, -74.0060),
    "Los Angeles": (34.0522, -118.2437),
    "Londyn": (51.5074, -0.1278),
    "Paryż": (48.8566, 2.3522),
    "Tokio": (35.6762, 139.6503),
    "Sydney": (-33.8688, 151.2093),
    "Kair": (30.0444, 31.2357),
    "Kapsztad": (-33.9249, 18.4241),
    "Bombaj": (19.0760, 72.8777),
    "Reykjavik": (64.1466, -21.9426),
    "Buenos Aires": (-34.6037, -58.3816),
    "Dubaj": (25.2048, 55.2708),
    "Nuuk": (64.1743, -51.7373)
}

#Etap I
class Wezel:
    def __init__(self, nazwa, pozycja_xyz=None, lat=None, lon=None):
        self.nazwa = nazwa
        self.pozycja_xyz = pozycja_xyz
        self.lat = lat
        self.lon = lon
        self.sasiedzi = []

    def dodaj_krawedz(self, sasiad, dystans):
        self.sasiedzi.append((sasiad, dystans))

class Graf:
    def __init__(self):
        self.wezly = []

    def dodaj_wezel(self, wezel):
        self.wezly.append(wezel)

def czy_ziemia_zaslania(p1, p2):
    R_ziemi = 6350.0
    d = p2 - p1
    a = np.dot(d, d)
    b = 2 * np.dot(p1, d)
    c = np.dot(p1, p1) - R_ziemi ** 2
    delta = b ** 2 - 4 * a * c

    if delta > 0:
        t1 = (-b - np.sqrt(delta)) / (2 * a)
        t2 = (-b + np.sqrt(delta)) / (2 * a)
        if (0 <= t1 <= 1) or (0 <= t2 <= 1):
            return True
    return False

#Etap II
def zbuduj_siec(limit_satelitow, nazwa_start, nazwa_cel):
    graf = Graf()
    ts = load.timescale()
    obecny_czas = ts.now()

    wszystkie_tle = load.tle_file('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle')
    satelity_tle = wszystkie_tle[:limit_satelitow]

    for sat in satelity_tle:
        geocentryczne = sat.at(obecny_czas)
        subpoint = geocentryczne.subpoint()

        obiekt_itrs = wgs84.latlon(subpoint.latitude.degrees, subpoint.longitude.degrees,
                                   elevation_m=subpoint.elevation.m)
        pozycja_xyz = obiekt_itrs.itrs_xyz.km

        nowy_wezel = Wezel(
            nazwa=sat.name,
            pozycja_xyz=pozycja_xyz,
            lat=subpoint.latitude.degrees,
            lon=subpoint.longitude.degrees
        )
        graf.dodaj_wezel(nowy_wezel)

    lat_s, lon_s = MIASTA[nazwa_start]
    lat_c, lon_c = MIASTA[nazwa_cel]

    pos_s = wgs84.latlon(lat_s, lon_s)
    pos_c = wgs84.latlon(lat_c, lon_c)

    wezel_start = Wezel(f"Połączenie z {nazwa_start}", pos_s.itrs_xyz.km, lat_s, lon_s)
    wezel_cel = Wezel(f"Połączenie do {nazwa_cel}", pos_c.itrs_xyz.km, lat_c, lon_c)

    graf.dodaj_wezel(wezel_start)
    graf.dodaj_wezel(wezel_cel)

    limit_ziemia_kosmos = 1200
    limit_kosmos_kosmos = 2500

    for i in range(len(graf.wezly)):
        for j in range(i + 1, len(graf.wezly)):
            w1 = graf.wezly[i]
            w2 = graf.wezly[j]
            dystans = np.linalg.norm(w1.pozycja_xyz - w2.pozycja_xyz)

            czy_polaczenie_z_ziemia = "Miasto" in w1.nazwa or "Miastio" in w2.nazwa
            aktywny_limit = limit_ziemia_kosmos if czy_polaczenie_z_ziemia else limit_kosmos_kosmos

            if dystans < aktywny_limit and not czy_ziemia_zaslania(w1.pozycja_xyz, w2.pozycja_xyz):
                w1.dodaj_krawedz(w2, dystans)
                w2.dodaj_krawedz(w1, dystans)

    return graf, wezel_start, wezel_cel

# Etap III
def znajdz_najkrotsza_trase(start, cel):
    dystanse = {start: 0}
    poprzednicy = {}
    kopiec = [(0, id(start), start)]
    odwiedzone = set()

    while kopiec:
        obecny_dystans, _, obecny_wezel = heapq.heappop(kopiec)

        if obecny_wezel == cel:
            break
        if obecny_wezel in odwiedzone:
            continue

        odwiedzone.add(obecny_wezel)

        for sasiad, waga_krawedzi in obecny_wezel.sasiedzi:
            if sasiad in odwiedzone:
                continue

            nowy_dystans = obecny_dystans + waga_krawedzi

            if sasiad not in dystanse or nowy_dystans < dystanse[sasiad]:
                dystanse[sasiad] = nowy_dystans
                poprzednicy[sasiad] = obecny_wezel
                heapq.heappush(kopiec, (nowy_dystans, id(sasiad), sasiad))

    trasa = []
    wezel_krok = cel

    if cel not in poprzednicy and cel != start:
        return None, float('inf')

    while wezel_krok:
        trasa.append(wezel_krok)
        wezel_krok = poprzednicy.get(wezel_krok)

    trasa.reverse()
    return trasa, dystanse.get(cel, float('inf'))

# Etap IV
def przeprowadz_testy_wydajnosci(nazwa_start, nazwa_cel):
    print("\nAnaliza")
    rozmiary_danych = [100, 300, 500, 800, 1200, 2000]
    czasy_wykonania_ms = []
    pingi_ms = []

    satelity_tle_pelne = load.tle_file('https://celestrak.org/NORAD/elements/gp.php?GROUP=starlink&FORMAT=tle')

    for N in rozmiary_danych:
        print(f"Obliczanie parametrów dla {N} satelitów")
        graf_testowy = Graf()

        lat_s, lon_s = MIASTA[nazwa_start]
        lat_c, lon_c = MIASTA[nazwa_cel]
        pos_s = wgs84.latlon(lat_s, lon_s)
        pos_c = wgs84.latlon(lat_c, lon_c)
        wezel_start = Wezel(f"Połączenie z {nazwa_start}", pos_s.itrs_xyz.km, lat_s, lon_s)
        wezel_meta = Wezel(f"Połączenie do {nazwa_cel}", pos_c.itrs_xyz.km, lat_c, lon_c)

        graf_testowy.dodaj_wezel(wezel_start)
        graf_testowy.dodaj_wezel(wezel_meta)

        ts = load.timescale()
        obecny_czas = ts.now()
        for sat in satelity_tle_pelne[:N]:
            pozycja_xyz = sat.at(obecny_czas).position.km
            graf_testowy.dodaj_wezel(Wezel(nazwa=sat.name, pozycja_xyz=pozycja_xyz))

        limit_ziemia_kosmos = 1200
        limit_kosmos_kosmos = 2500

        for i in range(len(graf_testowy.wezly)):
            for j in range(i + 1, len(graf_testowy.wezly)):
                w1 = graf_testowy.wezly[i]
                w2 = graf_testowy.wezly[j]
                dystans = np.linalg.norm(w1.pozycja_xyz - w2.pozycja_xyz)
                czy_polaczenie_z_ziemia = "Miasto" in w1.nazwa or "Miastio" in w2.nazwa
                aktywny_limit = limit_ziemia_kosmos if czy_polaczenie_z_ziemia else limit_kosmos_kosmos

                if dystans < aktywny_limit and not czy_ziemia_zaslania(w1.pozycja_xyz, w2.pozycja_xyz):
                    w1.dodaj_krawedz(w2, dystans)
                    w2.dodaj_krawedz(w1, dystans)

        start_time = time.perf_counter_ns()
        trasa, dystans = znajdz_najkrotsza_trase(wezel_start, wezel_meta)
        end_time = time.perf_counter_ns()

        czas_ms = (end_time - start_time) / 1_000_000
        czasy_wykonania_ms.append(czas_ms)

        if trasa:
            ping_rtt = (dystans / 300_000) * 1000 * 2
            pingi_ms.append(ping_rtt)
        else:
            pingi_ms.append(None)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle(f'Analiza struktur: {nazwa_start} -> {nazwa_cel}', fontsize=14)

    ax1.plot(rozmiary_danych, czasy_wykonania_ms, marker='o', color='red', linewidth=2)
    ax1.set_title('Czas wykonania algorytmu Dijkstry')
    ax1.set_xlabel('Rozmiar grafu (liczba satelitów)')
    ax1.set_ylabel('Czas obliczeń [ms]')
    ax1.grid(True, linestyle='--', alpha=0.7)

    valid_N = [N for N, p in zip(rozmiary_danych, pingi_ms) if p is not None]
    valid_ping = [p for p in pingi_ms if p is not None]

    ax2.plot(valid_N, valid_ping, marker='s', color='blue', linewidth=2)
    ax2.set_title('Optymalizacja ścieżki (opóźnienie RTT)')
    ax2.set_xlabel('Rozmiar grafu (liczba satelitów)')
    ax2.set_ylabel('Opóźnienie transmisji RTT [ms]')
    ax2.grid(True, linestyle='--', alpha=0.7)

    plt.tight_layout()
    plt.show()

# ETAP V
def rysuj_globus_3d(trasa):
    if not trasa:
        return

    lats = [w.lat for w in trasa]
    lons = [w.lon for w in trasa]
    names = [w.nazwa for w in trasa]
    kolory = ['blue'] + ['orange'] * (len(trasa) - 2) + ['green']

    fig = go.Figure()
    fig.add_trace(go.Scattergeo(
        lon=lons, lat=lats, mode='lines+markers+text', text=names, textposition="top center",
        textfont=dict(family="Arial", size=10, color="darkred"), line=dict(width=3, color='red'),
        marker=dict(size=8, color=kolory, line=dict(width=1, color='white')), name='Trasa LEO'
    ))

    fig.update_layout(
        title_text='Dynamiczny Routing LEO', title_x=0.5, showlegend=False,
        geo=dict(projection_type='orthographic', showcoastlines=True, coastlinecolor="black",
                 showland=True, landcolor="lightgray", showocean=True, oceancolor="aliceblue",
                 showcountries=True, countrycolor="darkgray"),
        margin={"r": 0, "t": 40, "l": 0, "b": 0}
    )
    fig.show()

# Etap VI
class AplikacjaSatelitarna:
    def __init__(self, root):
        self.root = root
        self.root.title("Panel sterowania")
        self.root.geometry("550x380")

        self.ostatnia_trasa = None
        lista_nazw = sorted(list(MIASTA.keys()))

        tk.Label(root, text="Projekt Routing Satelitarny LEO", font=("Arial", 16, "bold")).pack(pady=15)

        ramka_ustawien = tk.LabelFrame(root, text="Parametry", font=("Arial", 10, "bold"), padx=20,
                                       pady=10)
        ramka_ustawien.pack(fill="x", padx=20, pady=5)

        ramka_trasy = tk.Frame(ramka_ustawien)
        ramka_trasy.pack(pady=5)

        tk.Label(ramka_trasy, text="Skąd:", font=("Arial", 10)).grid(row=0, column=0, padx=5, sticky="e")
        self.combo_start = ttk.Combobox(ramka_trasy, values=lista_nazw, state="readonly", width=15, font=("Arial", 10))
        self.combo_start.set("Poznań")
        self.combo_start.grid(row=0, column=1, padx=5, pady=2)

        tk.Label(ramka_trasy, text="Dokąd:", font=("Arial", 10)).grid(row=0, column=2, padx=5, sticky="e")
        self.combo_cel = ttk.Combobox(ramka_trasy, values=lista_nazw, state="readonly", width=15, font=("Arial", 10))
        self.combo_cel.set("Nowy Jork")
        self.combo_cel.grid(row=0, column=3, padx=5, pady=2)

        ramka_input = tk.Frame(ramka_ustawien)
        ramka_input.pack(pady=5)
        tk.Label(ramka_input, text="Liczba satelitów:", font=("Arial", 10)).pack(side=tk.LEFT, padx=5)
        self.entry_satelitow = tk.Entry(ramka_input, font=("Arial", 10, "bold"), width=8, justify="center")
        self.entry_satelitow.insert(0, "2000")
        self.entry_satelitow.pack(side=tk.LEFT, padx=5)

        self.typ_raportu = tk.StringVar(value="prosty")
        ramka_radio = tk.Frame(ramka_ustawien)
        ramka_radio.pack(pady=5)
        tk.Radiobutton(ramka_radio, text="Raport podstawowy", variable=self.typ_raportu, value="prosty").pack(
            side=tk.LEFT, padx=15)
        tk.Radiobutton(ramka_radio, text="Raport zaawansowany", variable=self.typ_raportu, value="szczegolowy").pack(
            side=tk.LEFT, padx=15)

        ramka_przyciski = tk.Frame(root)
        ramka_przyciski.pack(pady=15)

        self.btn_start = tk.Button(ramka_przyciski, text="Przelicz Sieć", command=self.uruchom_w_tle, bg="#4CAF50",
                                   fg="white", font=("Arial", 11, "bold"), padx=5)
        self.btn_start.pack(side=tk.LEFT, padx=8)

        self.btn_mapa = tk.Button(ramka_przyciski, text="Wizualizacja 3D", command=self.wywolaj_mape, bg="#2196F3",
                                  fg="white", font=("Arial", 11, "bold"), state=tk.DISABLED, padx=5)
        self.btn_mapa.pack(side=tk.LEFT, padx=8)

        self.btn_wykres = tk.Button(ramka_przyciski, text="Wykresy", command=self.uruchom_testy_w_tle,
                                    bg="#FF9800", fg="white", font=("Arial", 11, "bold"), padx=5)
        self.btn_wykres.pack(side=tk.LEFT, padx=8)

        print("[SYSTEM] Panel Inicjalizacyjny gotowy. Dane wyjściowe przekierowane do terminala.")

    def wywolaj_mape(self):
        if self.ostatnia_trasa:
            rysuj_globus_3d(self.ostatnia_trasa)

    def uruchom_w_tle(self):
        self.btn_start.config(state=tk.DISABLED)
        self.btn_mapa.config(state=tk.DISABLED)
        watek = threading.Thread(target=self.wykonaj_obliczenia)
        watek.start()

    def wykonaj_obliczenia(self):
        try:
            limit = int(self.entry_satelitow.get())
            if limit < 20 or limit > 8000:
                raise ValueError("Można tylko podać liczbę od 20 do 8000.")

            miasto_s = self.combo_start.get()
            miasto_c = self.combo_cel.get()

            if miasto_s == miasto_c:
                messagebox.showwarning("Błąd", "Miasto początkowe i końcowe muszą być inne")
                return

            print(f"\nObliczanie ścieżki dla trasy: {miasto_s} -> {miasto_c} (Ilość satelitów = {limit})")

            t0 = time.perf_counter_ns()
            graf, start, meta = zbuduj_siec(limit, miasto_s, miasto_c)
            dt_budowy = (time.perf_counter_ns() - t0) / 1_000_000
            print(f"Macierz sąsiedztwa i graf stworzone w {dt_budowy:.1f} ms.")

            t1 = time.perf_counter_ns()
            trasa, dystans = znajdz_najkrotsza_trase(start, meta)
            dt_dijkstry = (time.perf_counter_ns() - t1) / 1_000_000

            if trasa:
                skoki = len(trasa) - 1
                rtt_ms = ((dystans / 300_000) * 1000) * 2
                self.ostatnia_trasa = trasa

                print(f"Trasa algorytmu Dijkstry. Łączny dystans: {dystans:.1f} km.")

                if self.typ_raportu.get() == "prosty":
                    tresc_raportu = (
                        f"Wynik raportu prostego\n\n"
                        f"Trasa: {miasto_s} ➔ {miasto_c}\n"
                        f"Rozmiar sieci: {limit} satelit\n"
                        f"Liczba przeskoków: {skoki}\n"
                        f"Dystans trasy: {dystans:.2f} km\n"
                        f"Szacowane opóźnienie (RTT): {rtt_ms:.2f} ms\n"
                    )
                else:
                    sprzet_upload_mbps = 25.0
                    tcp_limit_upload_mbps = (524288 / (rtt_ms / 1000)) / 1_000_000
                    realny_upload_mbps = min(tcp_limit_upload_mbps, sprzet_upload_mbps)

                    sprzet_download_mbps = 250.0
                    tcp_limit_download_mbps = (16777216 / (rtt_ms / 1000)) / 1_000_000
                    realny_download_mbps = min(tcp_limit_download_mbps, sprzet_download_mbps)

                    tresc_raportu = (
                        f"Wynik raportu zaawansowanego\n\n"
                        f"Parametry fizyczne:\n"
                        f"Relacja: {miasto_s} do {miasto_c}\n"
                        f"Przeskoki satelit: {skoki} | Dystans linku: {dystans:.2f} km\n"
                        f"Opóźnienie pętli RTT: {rtt_ms:.2f} ms\n\n"
                        f"Parametry wysyłania:\n"
                        f"Przepustowość sprzętowa: {sprzet_upload_mbps} Mbps\n"
                        f"Limit protokołu TCP: {tcp_limit_upload_mbps:.2f} Mbps\n"
                        f"Efektywny transfer danych: {realny_upload_mbps:.2f} Mbps\n\n"
                        f"Parametry pobierania:\n"
                        f"Przepustowość sprzętowa stacji: {sprzet_download_mbps} Mbps\n"
                        f"Limit dla TCP (2MB): {tcp_limit_download_mbps:.2f} Mbps\n"
                        f"Efektywny transfer danych: {realny_download_mbps:.2f} Mbps\n"
                    )

                messagebox.showinfo("Wynik", tresc_raportu)
                self.btn_mapa.config(state=tk.NORMAL)
            else:
                self.ostatnia_trasa = None
                print("Brak ciągłości geometrycznej.")
                messagebox.showwarning("Błąd",
                                       "Brak ciągłej trasy satelitarnej między wybranymi miastami dla danej liczby satelitów")

        except ValueError as e:
            messagebox.showerror("Błąd", f"Podaj poprawną liczbę satelitów:\n{e}")
        finally:
            self.btn_start.config(state=tk.NORMAL)

    def uruchom_testy_w_tle(self):
        self.btn_wykres.config(state=tk.DISABLED)
        watek = threading.Thread(target=self.wykonaj_testy)
        watek.start()

    def wykonaj_testy(self):
        try:
            miasto_s = self.combo_start.get()
            miasto_c = self.combo_cel.get()
            if miasto_s == miasto_c:
                messagebox.showwarning("Błąd", "Wybierz dwa różne miasta do wygenerowania wykresów")
                return
            przeprowadz_testy_wydajnosci(miasto_s, miasto_c)
        except Exception as e:
            print(f"Błedy z testów: {str(e)}")
        finally:
            self.btn_wykres.config(state=tk.NORMAL)


if __name__ == "__main__":
    okno = tk.Tk()
    app = AplikacjaSatelitarna(okno)
    okno.mainloop()
