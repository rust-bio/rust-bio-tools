const unpacker = new jsonm.Unpacker();
spec = unpacker.unpack(spec);
let changed_records = []
let delete_records = []
spec.datasets.main.forEach(function(record, j) {
    let variants = record.variants.split("/");
    if (variants.length > 1) {
        for (var i = 0; i < variants.length; i++) {
            let field = "variants" + (i+1);
            let layer_record = JSON.parse(JSON.stringify(record));
            layer_record[field] = variants[i];
            layer_record.count_variants = variants[i];
            changed_records.push(layer_record);
        }
        delete_records.push(j);
    } else {
        record.variants1 = record.variants;
        record.count_variants = record.variants;
    }
});
for (let i = 0; i < delete_records.length; i++) {
    spec.datasets.main.splice(delete_records[i]-i, 1);
}
changed_records.forEach(function(record) {
    spec.datasets.main.push(record);
});
let page_width =  $(window).width();
let matrix_width = Math.min(page_width - 740, samples*20);
if (matrix_width < 20 && samples > 2) {
    matrix_width = 40;
}
// add sort attribute to data
spec.vconcat[1].hconcat.forEach(function(ele) {
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

spec.vconcat[0].width = matrix_width;
vegaEmbed('#oncoprint', spec).then(function(result) {
    result.view.addEventListener('click', function(event, item) {
        if (item.datum.sample !== undefined) {
            window.location.href = '../details/' + item.datum.sample + '/' + gene + '.html';
        }
    });
});

window.addEventListener('resize', function(event){
    let page_width =  $(window).width();
    let matrix_width = Math.min(page_width - 740, samples*20);
    if (matrix_width < 20 && samples > 2) {
        matrix_width = 40;
    }
    spec.vconcat[0].width = matrix_width;
    spec.vconcat[1].hconcat.forEach(function(ele) {
        if (ele.layer !== undefined) {
            ele.layer.forEach(function(layer) {
                layer.width = matrix_width;
            })
        }
    });

    vegaEmbed('#oncoprint', spec).then(function(result) {
        result.view.addEventListener('click', function(event, item) {
            if (item.datum.sample !== undefined) {
                window.location.href = '../details/' + item.datum.sample + '/' + gene + '.html';
            }
        });
    });
});