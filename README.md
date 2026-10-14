# Symulator Routingu Satelitarnego LEO (dla konstelacji Starlink)

Aplikacja do analizy i wizualizacji 3D dynamicznego routingu w sieciach satelitów na niskiej orbicie Ziemi. Projekt symuluje przesyłanie danych między ziemskimi stacjami bazowymi z wykorzystaniem satelitów konstelacji Starlink, wyszukując najkrótszą ścieżkę komunikacyjną.

Ważnym elementem, jest to, że aplikacja jest projektem na studia, z takiego powodu zostały użyte inne algorytmy niż te, które powinny zostać użyte (np. K-Tree). Natomiast całe algorytmy były tworzone i analizowane pod wodzą **prowadzącego** zajęcia, przez co zaprogramowane algorytmy są w **100%** działające. 

## Główne funkcjonalności aplikacji

* **Pobieranie danych w czasie rzeczywistym**
  
  Aplikacja pobiera aktualne dane TLE satelitów Starlink z bazy Celestrak przy użyciu biblioteki `skyfield`.
* **Budowa grafu widoczności**
  
  Dynamiczne tworzenie węzłów i krawędzi z uwzględnieniem krzywizny Ziemi oraz maksymalnych zasięgów (1200 km dla połączeń Ziemia-Kosmos i 2500 km dla Kosmos-Kosmos).
* **Optymalizacja trasy**
  
  Implementacja algorytmu Dijkstry z wykorzystaniem kopca (`heapq`) do znalezienia najkrótszej ścieżki (najmniejszego opóźnienia) między miastami.
* **Szacowanie parametrów łącza**
  
  Obliczanie opóźnienia RTT oraz estymacja realnej przepustowości TCP dla pobierania i wysyłania danych.
* **Wizualizacja 3D**
  
  Generowanie interaktywnej Ziemi z wyznaczoną trasą przy użyciu biblioteki `plotly`.
* **Analiza wydajności**
  
  Moduł testujący czas wykonania algorytmu oraz opóźnienia w zależności od wielkości sieci (ilości satelitów), generujący wykresy w `matplotlib`.

## Użyte biblioteki

* **Język:** Python
* **GUI:** Tkinter
* **Obliczenia przestrzenne i astronomiczne:** `skyfield`, `numpy`
* **Wizualizacja:** `plotly`, `matplotlib`
* **Algorytmy:** `heapq`

## Użycie i instalacja

* Dla otwarcia pobranego folderu zalecam użycie środowiska **PyCharm**.

* Po uruchomieniu panelu wybierz z list rozwijanych miasto początkowe oraz miasto docelowe (np. Poznań -> Nowy Jork)
* Określ liczbę satelitów (od 20 do 8000), które mają zostać uwzględnione w routingu
* Wybierz typ raportu (podstawowy lub zaawansowany)
* Kliknij **Przelicz Sieć**.
* Po udanym wyznaczeniu trasy odblokuje się przycisk **Wizualizacja 3D**, który otworzy interaktywny model Ziemi w przeglądarce internetowej
* Przycisk **Wykresy** pozwala na wygenerowanie analizy  algorytmu dla różnych rozmiarów grafu (od 100 do 2000 satelitów)
