const modal = document.querySelector('.modal');
const input = modal.querySelector('input');

document.querySelector('.modal-close').onclick = () => modal.style.display = 'none';

window.addEventListener('click', e => e.target === modal ? modal.style.display = 'none' : null);

const populate_modal = txt => {
    input.value = txt;
    input.select();
};

const show_modal = () => {
    input.value = null;
    modal.style.display = 'block';
};

const copy_lsv_modal = (n, url) => {
    n.querySelectorAll('.copy-lsv')
        .forEach(c => c.onclick = () => {
            show_modal();
            json_ajax(url)
                .then(lsv_data => {
                    populate_modal(JSON.stringify(lsv_data).replace(/"/g, '\\"'))
                });
        });
};

const copy_lsv_modal_dropdown = (n, url) => {
    n.querySelectorAll('.copy-lsv')
        .forEach(c => c.onchange = () => {
            const grp = c.options[c.selectedIndex].getAttribute('value');
            if (grp) {
                show_modal();
                send_ajax(url, {'group_name': grp}).then(lsv_data => {
                    if(lsv_data != ''){
                        populate_modal(JSON.stringify(lsv_data).replace(/"/g, '\\"'))
                    }else{
                        populate_modal('No LSV data for the selected group');
                    }
                })
            }
        })
};