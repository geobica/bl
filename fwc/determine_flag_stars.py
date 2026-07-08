import json
import os

def remove_character_ranges(input_filename, output_filename, ranges):
    with open(input_filename, "r", encoding="utf-8") as f:
        text = f.read()

    keep = [True] * len(text)

    for start, end in ranges:
        start = max(0, start)
        end = min(len(text) - 1, end)
        if start <= end:
            for i in range(start, end + 1):
                keep[i] = False

    result = "".join(ch for ch, k in zip(text, keep) if k)

    with open(output_filename, "w", encoding="utf-8") as f:
        f.write(result)

def main():
    with open("svg_star_ref.json", "r", encoding="utf-8") as f:
        data = json.load(f)

    for element_name, element_data in data.items():
        input_filename = element_data["file"]
        ranges = [thing["svg_loc"] for thing in element_data["stars"]]

        base, ext = os.path.splitext(f"svg/original/{input_filename}")
        output_filename = f"svg/starless/{input_filename}"

        remove_character_ranges(
            f"svg/original/{input_filename}",
            output_filename,
            ranges,
        )

if __name__ == "__main__":
    main()