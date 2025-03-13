def measure_overlap(seq1, seq2):
    overlaps = [0]
    for i in range(1, len(seq2)):
        if seq2[:i] == seq1[-i:]: overlaps.append(i)
    return max(overlaps)
def measure_frag(seq1, seq2):
    overlap = measure_overlap(seq1, seq2)
    mer_size = len(seq2)
    hop = len(seq2) - overlap
    return dict(length=mer_size, step=hop, overlap=overlap)
def measure_frags(frags, labels):
    results = {}
    previous_frag = ''
    for frag, label in zip(frags, labels):
        results[frag] = measure_frag(previous_frag, frag)
        if previous_frag == '':
            results[frag]['step'] = None
            results[frag]['overlap'] = None
        try: label = bool(int(label))
        except:
            label = str(label)
            if label.lower() == 'true': label = True
            elif label.lower() == 'false': label = False
            elif len(label) <=0 : label = False
            else: label = True
        results[frag]['label'] = label
        previous_frag = frag
    return results