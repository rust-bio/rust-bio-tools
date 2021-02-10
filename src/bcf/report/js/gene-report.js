let decompressed = LZString.decompressFromUTF16(spec);
const unpacker = new jsonm.Unpacker();
spec = unpacker.unpack(JSON.parse(decompressed));
let changed_records = []
let delete_records = []
spec.datasets.main.forEach(function(record, j) {
    let variant = record.variant;
    record.variants1 = [record.variant];
    record.count_variants = [record.variant];
});
for (let i = 0; i < delete_records.length; i++) {
    spec.datasets.main.splice(delete_records[i]-i, 1);
}
changed_records.forEach(function(record) {
    spec.datasets.main.push(record);
});
let page_width =  $(window).width();
let matrix_width = Math.min(page_width - 740, samples*20);
if (matrix_width < 20 && samples >= 2) {
    matrix_width = 40;
} else if (samples == 1) {
    matrix_width = 20;
}

let pl = spec.vconcat.length - 1;
// add sort attribute to data
spec.vconcat[pl].hconcat.forEach(function(ele) {
    if (ele.layer === undefined) {
        ele.encoding.y.sort = order;
    } else {
        ele.layer.forEach(function(layer) {
            layer.encoding.y.sort = order;
            layer.width = matrix_width;
        })
    }
})

var gene = $(document).attr('title');

for (let z = 0; z < pl; z++) {
    spec.vconcat[z].width = matrix_width;
}
vegaEmbed('#oncoprint', spec).then(function(result) {
    result.view.addEventListener('click', function(event, item) {
        if (item.datum.sample !== undefined  && item.datum.alteration !== undefined) {
            window.location.href = '../details/' + item.datum.sample + '/' + gene + '.html';
        }
    });
});

window.addEventListener('resize', function(event){
    let page_width =  $(window).width();
    let matrix_width = Math.min(page_width - 740, samples*20);
    if (matrix_width < 20 && samples >= 2) {
        matrix_width = 40;
    } else if (samples == 1) {
        matrix_width = 20;
    }
    for (let z = 0; z < pl; z++) {
        spec.vconcat[z].width = matrix_width;
    }
    spec.vconcat[pl].hconcat.forEach(function(ele) {
        if (ele.layer !== undefined) {
            ele.layer.forEach(function(layer) {
                layer.width = matrix_width;
            })
        }
    });

    vegaEmbed('#oncoprint', spec).then(function(result) {
        result.view.addEventListener('click', function(event, item) {
            if (item.datum.sample !== undefined  && item.datum.alteration !== undefined) {
                window.location.href = '../details/' + item.datum.sample + '/' + gene + '.html';
            }
        });
    });
});
